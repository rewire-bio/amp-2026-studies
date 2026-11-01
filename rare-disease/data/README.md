# Rare-disease data

Every input and intermediate file of the case benchmark. With these files committed, `case_benchmark.py` reruns
without downloading anything. All files are tab-separated with a header row, on GRCh38, with chr-prefixed
chromosomes, unless noted.

## Layout

```
data/                        SNV study (200 cases) and shared inputs
  indel/                     indel study (120 cases)
  cancer/                    inherited cancer study (140 cases), its own panel, regions and backgrounds
  replication/               replication study (500 cases), 50 new backgrounds
  kgp/                       1000 Genomes genotypes and lookup region files for the SNV and indel studies
  spliceai/                  SpliceAI run on every indel
  logs/                      run logs
```

Each study folder has the same core files:

| File | Made by | Contents |
|---|---|---|
| `cases.tsv` | `build_cases.py` | One row per variant per case: `case`, `sample` (1000 Genomes person), `chrom`, `pos`, `ref`, `alt`, `diagnostic` (1 for the inserted ClinVar variant), `group` (diagnostic consequence group), `gene`, `clinvar_id` |
| `case_variants.tsv` | `build_cases.py` | The distinct variants across all cases (`chrom`, `pos`, `ref`, `alt`): the list every score is looked up for |
| `scores_<source>.tsv` | `lookup_scores.py` | One score per variant, keyed by `chrom`, `pos`, `ref`, `alt`; a variant missing from the file has no score |
| `consequence_indels.tsv` | `annotate_indels.py` | `indel_class` of every indel from MANE Select (indel and cancer studies) |

Score files and their columns:

| File | Columns | Notes |
|---|---|---|
| `scores_cadd.tsv` | `cadd_phred`, `cadd_consequence`, `gene`, `spliceai_max` | CADD v1.7 inclAnno; SpliceAI is the largest of its four delta scores |
| `scores_revel.tsv` | `revel` | REVEL v1.3, highest over transcripts; missense only |
| `scores_am.tsv` | `alphamissense` | AlphaMissense, highest over transcripts; missense only |
| `scores_gpnmsa.tsv` | `gpn_msa` | GPN-MSA, sign flipped so higher means more damaging |
| `scores_gpnstar_v.tsv`, `scores_gpnstar_m.tsv` | `gpnstar_v`, `gpnstar_m` | GPN-Star vertebrate and mammal, calibrated, sign flipped |
| `scores_avi.tsv` | `avi`, `avi_phred` | AlphaGenome AVI raw score and its PHRED-style scale; SNVs only |
| `scores_spliceai_run.tsv` | `spliceai_run` | SpliceAI 1.3.1 run ourselves on indels |

Evo 2 scores are not here; they are in [`../evo2/run/`](../evo2/run/) and [`../evo2/run3/`](../evo2/run3/).

## Shared inputs (top level)

| File | Made by | Contents |
|---|---|---|
| `panelapp_285.json` | `prepare_data.sh` | PanelApp panel 285, Intellectual disability, v11.45 |
| `MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` | `prepare_data.sh` | MANE Select transcripts |
| `regions.tsv`, `regions_nochr.tsv`, `region_genes.tsv` | `regions.py` | Exons of the panel's green genes plus 20 bp, merged, 1-based inclusive; the same without `chr`; the gene each region came from |
| `clinvar_plp_panel.vcf.gz`, `clinvar_indel_panel.tsv` | `prepare_data.sh` | ClinVar records in the regions (fileDate 2026-10-11) |
| `clinvar_plp2.tsv`, `clinvar_indel_plp2.tsv` | `clinvar_plp.py` | ClinVar P/LP with 2+ review stars, one row per allele, most severe consequence: the pools diagnostic variants are drawn from |
| `all_samples.txt`, `samples50.txt` | `prepare_data.sh` | All 3,202 1000 Genomes high-coverage samples; the 50 drawn as backgrounds (seed 20261016) |
| `gnomad_af.tsv` | `gnomad_af.py` | For every variant in every study: `gnomad_found`, `gnomad_af` and `gnomad_af_grpmax` (highest of afr, amr, eas, mid, nfe, sas), from gnomAD v3.1.2 genomes via MyVariant.info |
| `cases_run1.tsv`, `case_variants_run1.tsv` | `build_cases.py` | An earlier build of the 200 SNV cases: the same 50 people, different diagnostic variants. Evo 2 run 1 scored these variants, so they are kept to trace `evo2/run/scores.tsv`; the extra variants of the final build are in `evo2/run/scores_extra.tsv` |

## kgp/

| File | Contents |
|---|---|
| `chrN.bcf` | 1000 Genomes high-coverage genotypes of the 50 background people in the regions (`fetch_kgp.sh`) |
| `regions_chrN.tsv` | The regions on that chromosome, as passed to bcftools |
| `cadd_regions_*`, `gpn_msa_regions*`, `avi_regions.tsv` | Position lists for the remote and local tabix lookups (`lookup_scores.py`) |

`cancer/kgp/` and `replication/kgp/` hold the same genotype files for those studies' backgrounds. The cancer
genotypes cover 21 chromosomes: panel 245 has no gene on chr8 or chr21.

## cancer/ and replication/

`cancer/` also has the cancer panel (`panelapp_245.json`), its regions and ClinVar pools (`clinvar_all_plp2.tsv`
combines SNVs and indels). `replication/` also has the 1000 Genomes pedigree and related-sample files used to pick 50
new, unrelated people (seed 20261023), and their `samples50.txt`.

## spliceai/

`indels.vcf` is every indel in the indel and cancer studies, as given to SpliceAI 1.3.1 (`run_spliceai.sh`), and
`indels.out.vcf` is the annotated output that `lookup_scores.py spliceai_run` reads.

## logs/

| File | Contents |
|---|---|
| `lookup-snv.log` | Score lookups for the 200 SNV cases: variants scored per source |
| `replication-fetch-kgp.log`, `replication-lookup-*.log` | Genotype fetch and score lookups for the replication |
| `benchmark-gnomad/<study>-<filter>.log` | Console output of each frequency-filtered `case_benchmark.py` run (`af001` is AF below 0.1%, `af01` below 1%, `grp001` grpmax below 0.1%) |

## Not committed

Downloaded by the scripts when needed, and listed in `.gitignore`:

- `avi/`: the AlphaGenome AVI all-SNV table, 82 GB unpacked (a symlink to an external disk is fine);
- `am/`: the AlphaMissense hg38 table, 640 MB;
- `clinvar.vcf.gz`: the full ClinVar VCF, 200 MB;
- `clinvar_snv_panel.tsv`: 167 MB, over GitHub's file size limit, rebuilt by `prepare_data.sh`;
- `*.tbi`, `*.csi`: tabix indexes of remote files.
