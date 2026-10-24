#!/usr/bin/env bash
# Download and prepare the public inputs for the rare-disease benchmark, in order. Run from rare-disease.
# Versions used for the talk: PanelApp 285 Intellectual disability v11.45, MANE v1.5, ClinVar fileDate 2026-10-11,
# 1000 Genomes high-coverage 3,202-sample phased panel (20220422).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data/kgp
cd data
KGP=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV

# 1. Panel, MANE exons, ClinVar.
curl -sSfLO https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz
curl -sSfLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -sSfLO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
curl -sSf "https://panelapp.genomicsengland.co.uk/api/v1/panels/285/" > panelapp_285.json

# 2. Target regions: MANE Select exons of green genes, plus 20 bp.
python3 ../regions.py
sed 's/^chr//' regions.tsv > regions_nochr.tsv

# 3. ClinVar pathogenic or likely pathogenic SNVs and indels (up to 50 bp) with at least two review stars, in the
#    regions. The review status is matched exactly in Python: bcftools treats commas in the string as OR.
FIELDS='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/MC\t%INFO/GENEINFO\t%INFO/CLNDN\n'
bcftools query -R regions_nochr.tsv -i 'INFO/CLNVC="single_nucleotide_variant"' -f "$FIELDS" clinvar.vcf.gz > clinvar_snv_panel.tsv
bcftools query -R regions_nochr.tsv -i 'INFO/CLNVC="Deletion" || INFO/CLNVC="Insertion" || INFO/CLNVC="Duplication" || INFO/CLNVC="Indel"' \
  -f "$FIELDS" clinvar.vcf.gz > clinvar_indel_panel.tsv
python3 ../clinvar_plp.py clinvar_snv_panel.tsv clinvar_plp2.tsv
python3 ../clinvar_plp.py clinvar_indel_panel.tsv clinvar_indel_plp2.tsv

# 4. Fifty random 1000 Genomes individuals, then their genotypes in the regions.
bcftools query -l "$KGP/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" > all_samples.txt
python3 -c "
import random; s=open('all_samples.txt').read().split(); random.Random(20261016).shuffle(s)
open('samples50.txt','w').write('\n'.join(sorted(s[:50]))+'\n')"
printf 'chr%s\n' {1..22} X | xargs -P 8 -n 1 ../fetch_kgp.sh

# 5. AlphaMissense GRCh38 (CC BY-NC-SA 4.0; kept local).
mkdir -p am && curl -sSfL -o am/AlphaMissense_hg38.tsv.gz https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz

# Then: python3 ../build_cases.py && python3 ../lookup_scores.py && python3 ../evo2/make_windows.py

# Replication (REPLICATION.md): 50 new backgrounds from the 2,504 unrelated samples, leaving out the first 50 and
# their families, then 500 new cases scored the same way.
mkdir -p replication && cd replication
curl -sSfLO "${KGP%/working/*}/20130606_g1k_3202_samples_ped_population.txt"
curl -sSfLO "${KGP%/working/*}/1000G_698_related_high_coverage.sequence.index"
cd ..
python3 -c "
import random
ped = [l.split() for l in open('replication/20130606_g1k_3202_samples_ped_population.txt')][1:]
fam = {r[1]: r[0] for r in ped}
related = {l.split('\t')[9] for l in open('replication/1000G_698_related_high_coverage.sequence.index') if not l.startswith('#')}
old = open('samples50.txt').read().split()
s = [x for x in open('all_samples.txt').read().split() if x not in related and fam[x] not in {fam[o] for o in old}]
random.Random(20261023).shuffle(s)
open('replication/samples50.txt','w').write('\n'.join(sorted(s[:50]))+'\n')"
printf 'chr%s\n' {1..22} X | xargs -P 8 -I{} ../fetch_kgp.sh {} replication/samples50.txt replication/kgp
# Then: python3 ../build_cases.py replication && python3 ../lookup_scores.py --study replication --reuse \
#   && python3 ../case_benchmark.py --study replication
