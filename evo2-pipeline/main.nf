#!/usr/bin/env nextflow
// Add an Evo 2 7B zero-shot score to a case's candidate variants, next to existing predictors, and rank them.
//
//   nextflow run main.nf -profile replay    # laptop: Evo 2 scores replayed from a recorded GPU run
//   nextflow run main.nf -profile gpu       # GPU host: scores computed with the pinned Evo 2 container
//
// Steps: normalise, attach precomputed scores (CADD, SpliceAI), build 8,192 bp windows for SNVs,
// score them (or replay recorded scores), attach the Evo 2 score, rank the case.
nextflow.enable.dsl = 2

params.case_vcf      = "${projectDir}/test/case_brca1.vcf"
params.reference     = "${projectDir}/resources/GRCh37_chr17.fa.gz"
params.cadd          = "${projectDir}/resources/cadd_v1.7.tsv.gz"
params.spliceai      = "${projectDir}/resources/spliceai_max.tsv.gz"
params.evo2_recorded = null      // replay mode: a recorded scores table (id, ref_score, var_score, delta)
params.window        = 8192
params.outdir        = "results"

process NORMALISE {
    input:
    path vcf
    path ref
    path ref_index

    output:
    path "normalised.vcf.gz"

    script:
    """
    bcftools norm --fasta-ref ${ref} --multiallelics -any --check-ref x ${vcf} -Oz -o normalised.vcf.gz
    """
}

process ANNOTATE_PRECOMPUTED {
    input:
    path vcf
    tuple path(cadd), path(cadd_index)
    tuple path(spliceai), path(spliceai_index)

    output:
    path "precomputed.vcf.gz"

    script:
    """
    cat > header.txt <<'EOF'
    ##INFO=<ID=CADD_RAW,Number=1,Type=Float,Description="CADD v1.7 raw score (GRCh37)">
    ##INFO=<ID=SPLICEAI_MAX,Number=1,Type=Float,Description="Largest SpliceAI masked delta score">
    EOF
    sed -i.bak 's/^ *//' header.txt
    bcftools annotate -a ${cadd} -c CHROM,POS,REF,ALT,CADD_RAW -h header.txt ${vcf} -Oz -o with_cadd.vcf.gz
    bcftools annotate -a ${spliceai} -c CHROM,POS,REF,ALT,SPLICEAI_MAX with_cadd.vcf.gz -Oz -o precomputed.vcf.gz
    """
}

process MAKE_WINDOWS {
    input:
    path vcf
    path ref
    path ref_index

    output:
    path "windows.tsv"

    script:
    """
    bcftools view --types snps ${vcf} | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' > snvs.tsv
    make_windows.py --snvs snvs.tsv --reference ${ref} --window ${params.window} --out windows.tsv
    """
}

process SCORE_EVO2 {
    label 'gpu'

    input:
    path windows

    output:
    path "evo2_scores.tsv"

    script:
    """
    score_evo2.py --windows ${windows} --out evo2_scores.tsv
    """
}

process REPLAY_EVO2 {
    input:
    path windows
    path recorded

    output:
    path "evo2_scores.tsv"

    script:
    """
    replay_evo2.py --windows ${windows} --recorded ${recorded} --out evo2_scores.tsv
    """
}

process ANNOTATE_EVO2 {
    publishDir params.outdir, mode: 'copy'

    input:
    path vcf
    path scores

    output:
    path "annotated.vcf.gz"

    script:
    """
    printf '##INFO=<ID=EVO2_LOF,Number=1,Type=Float,Description="Evo 2 7B zero-shot score: reference minus alternate mean log-likelihood per base, ${params.window} bp window. Higher is more damaging.">\\n' > header.txt
    tail -n +2 ${scores} | cut -f1-4,7 | sort -k1,1 -k2,2n | bgzip > evo2.tsv.gz
    tabix -s1 -b2 -e2 evo2.tsv.gz
    bcftools annotate -a evo2.tsv.gz -c CHROM,POS,REF,ALT,EVO2_LOF -h header.txt ${vcf} -Oz -o annotated.vcf.gz
    """
}

process RANK_CASE {
    publishDir params.outdir, mode: 'copy'

    input:
    path vcf

    output:
    path "ranked.tsv"

    script:
    """
    bcftools query -H -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CADD_RAW\\t%INFO/SPLICEAI_MAX\\t%INFO/EVO2_LOF\\n' ${vcf} > scores.tsv
    rank_case.py --scores scores.tsv --out ranked.tsv
    """
}

workflow {
    ref       = file(params.reference)
    ref_index = files("${params.reference}.{fai,gzi}")
    normalised  = NORMALISE(file(params.case_vcf), ref, ref_index)
    precomputed = ANNOTATE_PRECOMPUTED(normalised,
                                       tuple(file(params.cadd), file("${params.cadd}.tbi")),
                                       tuple(file(params.spliceai), file("${params.spliceai}.tbi")))
    windows = MAKE_WINDOWS(precomputed, ref, ref_index)
    scores  = params.evo2_recorded ? REPLAY_EVO2(windows, file(params.evo2_recorded)) : SCORE_EVO2(windows)
    RANK_CASE(ANNOTATE_EVO2(precomputed, scores))
}
