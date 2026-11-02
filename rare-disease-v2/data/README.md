# Data

Every file the second run reads or writes under `data/`. Tables are tab-separated with a header, GRCh38, with
chr-prefixed chromosomes. The variant key is `chrom, pos, ref, alt` (VCF positions). Scripts run from
`rare-disease-v2/`; the main README lists the order.

## Committed

### People and genotypes

| File | Made by | Columns and contents |
|---|---|---|
| `people.tsv` | `people.py` | One row per person in both sets. `sample` (1000 Genomes ID); `set` (`main`: the first run's 50; `rep`: the 50 new replication people); `sex` (`male` or `female`, from the 1000 Genomes pedigree); `gnomad_id` (the ID in gnomAD's HGDP + 1KG metadata, a trailing A or B stripped when matching); `gnomad_release` (`true` if gnomAD released the sample); `gnomad_pop` (gnomAD's inferred group); `family` (gnomAD IDs of the released family members, comma-separated, the person included if released) |
| `rep/samples50.txt` | `people.py` | The 50 replication people, one ID per line, sorted |
| `rep/kgp/chrN.bcf` | `fetch_rep_kgp.sh` (runs `../rare-disease/fetch_kgp.sh`) | 1000 Genomes high-coverage phased genotypes (20220422) of the 50 replication people in the panel regions. The main set's genotypes are read from `../rare-disease/data/kgp/` |
| `rep/kgp/regions_chrN.tsv` | `fetch_rep_kgp.sh` | The regions fetched on that chromosome (from `../rare-disease/data/regions.tsv`): `chrom`, `start`, `end`, no header |
| `main/person_variants.tsv`, `rep/person_variants.tsv` | `backgrounds.py` | Each person's variants with 1000 Genomes AF below 1% in the regions: `sample`, key, `zygosity` (`het`, `hom`, or `hemi` for a haploid call on chrX outside the PARs), `kgp_af` (1000 Genomes AF) |

### Frequencies

| File | Made by | Columns and contents |
|---|---|---|
| `gnomad_full.tsv` | `gnomad_af.py full` | Full gnomAD v3.1.2 genome counts from MyVariant.info for every person variant and the first run's diagnostics: key; `found` (1 if MyVariant.info returned a gnomAD genome record, which can have AC 0; `case_benchmark.py` counts a diagnostic as found only if `ac` is above 0); `ac`, `an` (overall); `ac_<g>`, `an_<g>` for g in afr, amr, eas, nfe, sas; `af` (overall AF); `grpmax` (highest AF over those five groups) |
| `main/person_variant_af.tsv`, `rep/person_variant_af.tsv` | `gnomad_af.py af` | One row per person variant: `sample`, key; `af_fam`, `grpmax_fam` (AF and grpmax with the person and their released family left out); `af_no1kg`, `grpmax_no1kg` (all released 1000 Genomes samples left out); `in_subset` (1 if the HGDP + 1KG subset file has a row at the variant) |
| `pool/gnomad_full.tsv` | `pool.py` | The same layout as `gnomad_full.tsv`, for the ClinVar candidate pools; `found` is `True` or `False` |

### Inheritance and candidate pools

| File | Made by | Columns and contents |
|---|---|---|
| `gene_moi.tsv` | `moi.py` | One row per MANE gene in the panel regions: `gene` (MANE symbol); `ensg`; `moi` (`mono`, `biallelic`, `xl_mono`, `xl_bi` or `unknown`); `panelapp_moi` (PanelApp's mode of inheritance string); `panelapp_gene` and `match` (`ensg`, `symbol`, `alias` or `none`), for checking |
| `pool/snv_pool.tsv` | `pool.py` | ClinVar P/LP 2+ star SNVs in the regions: ClinVar's `chrom`, `pos`, `clinvar_id`, `ref`, `alt`, `clnsig`, `revstat`, `consequence`, `clinvar_gene` (ClinVar's gene), `disease`; then `class` (section 4 class from CADD), `group` (diagnostic group or `none`), `cadd_consequence`, `cadd_mane` (1 if CADD has a MANE row), `grpmax` (full gnomAD), `found`, `gene` (MANE gene of the region), `moi`, `eligible` (1 if grpmax is below 0.1%) |
| `pool/indel_pool.tsv` | `pool.py` | The same for indels up to 50 bp, without `cadd_consequence` and `cadd_mane`; `class` from `annotate_indels.py`. Indels with REF and ALT both longer than 1 bp are left out |
| `pool/consequence_indels.tsv` | `pool.py` (calls `annotate_indels.py`) | key, `indel_class` (`splice`, `frameshift`, `inframe`, `utr`, `intronic` or `other`) for the indel pool |

### Cases

Three studies: `snv` (200 SNV cases from the main set), `indel` (120 indel cases from the main set) and
`replication` (500 SNV cases from the rep set).

| File | Made by | Columns and contents |
|---|---|---|
| `<study>/cases.tsv` | `build_cases.py` | One row per case variant. `case`, `sample`, key; `diagnostic` (1 for the diagnostic, 0 for background); `group` and `class` (diagnostic only; background class is assigned from CADD at analysis time); `gene` (every MANE region gene at the variant, comma-separated; for a diagnostic its pool gene first); `moi` (each gene's mode, same order); `zygosity`; `af`, `grpmax` (background: family left out; diagnostic: full gnomAD); `grpmax_no1kg` (background: all 1000 Genomes left out; diagnostic: equal to `grpmax`); `clinvar_id` (diagnostic only). The background holds SNVs and indels in every study, because the inheritance filter counts both; only the study's own type is ranked |
| `<study>/case_variants.tsv` | `build_cases.py` | The distinct keys across the study's cases |
| `snv/consequence_indels.tsv`, `indel/consequence_indels.tsv` | `annotate_indels.py --study <study>` | key, `indel_class`, for the indels in the study's cases (in the snv study, the background indels) |

### Scores

One file per source per study, written by `lookup_scores.py --study <study>`, one row per key found. A key with no row
has no score, and a missing score ranks last.

| File | Columns | Source |
|---|---|---|
| `<study>/scores_am.tsv` | key, `alphamissense` | AlphaMissense hg38 table, max over transcripts; SNVs |
| `<study>/scores_avi.tsv` | key, `avi`, `avi_phred` | AlphaGenome AVI all-SNV table (raw score and PHRED-scaled); SNVs |
| `<study>/scores_gpnmsa.tsv` | key, `gpn_msa` | GPN-MSA hg38 scores, -LLR; SNVs |
| `<study>/scores_gpnstar_v.tsv`, `scores_gpnstar_m.tsv` | key, `gpnstar_v` or `gpnstar_m` | GPN-Star vertebrate (v100) and mammal (m447), -llr_calibrated; SNVs |
| `<study>/scores_evee.tsv` | key, `evee` (pathogenicity), `evee_label` (EVEE's ClinVar label), `evee_variation_id` | EVEE (Zenodo 19701997, CC BY 4.0), ClinVar variants only, joined at EVEE pos equal to VCF POS minus 1 |
| `indel/scores_spliceai_run.tsv` | key, `spliceai_run` (max of the four SpliceAI deltas) | SpliceAI 1.3.1, run by `run_spliceai.sh`; indels |
| `spliceai/indels.out.vcf` | VCF with SpliceAI's `INFO/SpliceAI` field; `INFO` is `.` where SpliceAI gives no score | `run_spliceai.sh`, on the indels of `indel/case_variants.tsv` |

AlphaMissense, AVI and the GPN scores do not depend on annotation, so rows are copied from the first run where it
looked the key up. The SNV studies' SpliceAI score is `spliceai_max` from CADD's annotation (below).

## Not committed

| File | Why | How to get it back |
|---|---|---|
| `<study>/scores_cadd.tsv`, `pool/scores_cadd.tsv` | CADD is free for non-commercial use but not licensed for redistribution | `lookup_scores.py ... cadd` (remote tabix, no download); `pool.py` for the pool. Columns: key, `cadd_phred`, `cadd_consequence` (most severe over MANE Select rows, or all rows if none), `cadd_consequence_any` (most severe over all rows), `gene`, `mane` (1 if a MANE row exists), `spliceai_max` (largest of CADD's four SpliceAI deltas) |
| `<study>/scores_revel.tsv` | REVEL is free for non-commercial use but not licensed for redistribution | `lookup_scores.py ... revel`, from the REVEL v1.3 zip in `../../evo2-cloud-operations/data/raw/current-predictors/`. Columns: key, `revel` (max over transcripts); missense SNVs |
| `evee/evee_slim.parquet` | A 73 MB slim copy (ten columns) of the 35 GB EVEE table; kept local to keep the repository small. The per-case EVEE scores are committed | Built by `lookup_scores.py ... evee` on first use, by byte-range reads of the Zenodo shards; `results/evee_source.md` records the record version and checksums |
| `gnomad_subset/` | A local cache of downloads and intermediate files: the HGDP + 1KG sample metadata, the 1000 Genomes pedigree and sample lists, tabix indexes, slim VCFs of gnomAD's genotype calls at the person-variant positions (`dense/`), per-row adj counts (`counts_chrN.tsv`), fetch logs | `people.py` and `gnomad_af.py fetch` and `counts` (the fetch takes about an hour) |
| `spliceai/hg38.fa`, `hg38.fa.fai` | 3 GB, and already on disk in a sibling repository | `gunzip -c ../../evo2-size-by-use-case/data/raw/hg38.fa.gz > data/spliceai/hg38.fa`, then `samtools faidx` |
| `*.tbi`, `*.csi` | Index files | Rebuilt by the scripts that need them |
