# Work plan

Tick each task when its output exists. Design: [`PROTOCOL.md`](PROTOCOL.md). Why: [`PROPOSAL.md`](PROPOSAL.md).

Code is copied from [`../rare-disease/`](../rare-disease/) at tag `amp-2026-v4` and changed here; each script's
docstring names the file it came from. `../rare-disease/` code and data are not edited (only its README.md and
REPLICATION.md get notes pointing here). Large inputs (1000 Genomes genotypes for the first run's people, ClinVar,
MANE, PanelApp, the AVI table, AlphaMissense) are read from `../rare-disease/data/`. Per-variant scores that do not
depend on annotation (REVEL, AlphaMissense, AVI, GPN-MSA, GPN-Star) are reused from there; CADD is looked up again
because the consequence rule changed.

Run everything from `rare-disease-v2/` with `uv run --with-requirements requirements.txt python <script>`.

## Files

All tables are tab-separated with a header, GRCh38, chr-prefixed. The variant key is `chrom, pos, ref, alt`.

| File | Made by | Contents |
|---|---|---|
| `data/people.tsv` | `people.py` | `sample`, `set` (main: the first run's 50; rep: 50 new), `sex` (male or female), `gnomad_id`, `gnomad_release`, `gnomad_pop`, `family` (released family members' gnomAD IDs, comma-separated, the person included if released) |
| `data/rep/kgp/*.bcf` | `fetch_kgp.sh` (first run's) | genotypes of the new people in the panel regions |
| `data/<set>/person_variants.tsv` | `backgrounds.py` | `sample`, key, `zygosity` (het, hom or hemi), `kgp_af`: each person's variants with 1000 Genomes AF below 1% in the regions |
| `data/gnomad_full.tsv` | `gnomad_af.py` | key, `found`, `ac`, `an`, `ac_<g>`, `an_<g>` for g in afr, amr, eas, nfe, sas, `af`, `grpmax` (full gnomAD v3.1.2) for every person variant and pool variant |
| `data/<set>/person_variant_af.tsv` | `gnomad_af.py` | `sample`, key, `af_fam`, `grpmax_fam` (family left out), `af_no1kg`, `grpmax_no1kg` (all 1000 Genomes left out), `in_subset` (row in the HGDP + 1KG file) |
| `data/gene_moi.tsv` | `moi.py` | `gene` (MANE symbol), `ensg`, `moi` (mono, biallelic, xl_mono, xl_bi or unknown), `panelapp_moi` |
| `data/pool/snv_pool.tsv`, `data/pool/indel_pool.tsv` | `pool.py` | the ClinVar P/LP 2+ star candidates with `class` (protocol section 4), `grpmax` (full gnomAD), `gene`, `moi` |
| `data/<study>/cases.tsv` | `build_cases.py` | `case`, `sample`, key, `diagnostic`, `group`, `class`, `gene`, `moi`, `zygosity`, `af`, `grpmax`, `grpmax_no1kg`, `clinvar_id`; study is snv, indel or replication |
| `data/<study>/case_variants.tsv` | `build_cases.py` | distinct keys across the study's cases |
| `data/<study>/scores_<source>.tsv` | `lookup_scores.py` | as in the first run, plus `scores_evee.tsv` (`evee`, `evee_label`) |
| `data/<study>/consequence_indels.tsv` | `annotate_indels.py` | `indel_class` |
| `results/*.json` | `case_benchmark.py` | results |

## Phase 0: set-up

- [x] `PROPOSAL.md`, `PROTOCOL.md`, `PLAN.md`, `README.md`; links from the top-level README, GUIDE and the first run
- [x] Design review: feasibility probes (gnomAD, PP3 thresholds, EVEE, ClinVar pool) and three critics; corrections
  folded into the protocol (deviation log, section 13)

## Phase 1: people and frequencies (issue #1)

- [x] `people.py`: sex and family from the 1000 Genomes pedigree; gnomAD release status and group from the HGDP +
  1KG sample metadata; 50 new replication people (protocol section 5).
- [x] Genotypes of the new people: `fetch_rep_kgp.sh` (runs `../rare-disease/fetch_kgp.sh` with retries) into `data/rep/kgp/`.
- [x] `backgrounds.py`: `person_variants.tsv` for both sets.
- [x] `gnomad_af.py`: full counts from MyVariant.info; family and all-1000-Genomes counts from the HGDP + 1KG file;
  the checks in protocol section 1.
- [x] `results/gnomad_overlap.md`: the issue #1 table redone with the new frequencies, and the median background
  per case after the grpmax filter with the full, family-left-out and all-1000-Genomes-left-out frequencies.

## Phase 2: inheritance, pool and code

- [x] `moi.py`: `gene_moi.tsv`.
- [x] `annotate_indels.py`: takes `--variants` and `--out`, so it annotates the indel pool and the case variants.
- [ ] `lookup_scores.py`: cadd (MANE consequence rule), revel, am, avi, spliceai_run, gpnmsa, gpnstar, evee; reuses
  first-run rows; scores the pool for CADD.
- [x] `pool.py`: classified, frequency-checked, inheritance-tagged pools. Count positions per class against the
  quotas (200 and 500 together for SNVs, excluding the first run's diagnostic positions for the 500).
- [x] `build_cases.py`: protocol sections 2, 3 and 5.
- [ ] `case_benchmark.py`: protocol sections 7 to 11.
- [x] `run_spliceai.sh`, `evo2/make_windows.py` (`--fasta` required, `--indels-only`).
- [ ] Test: run `case_benchmark.py` with and without one secondary score; every row not involving it must match.

## Phase 3: build and first results

- [x] hg38 for SpliceAI: decompress `../../evo2-size-by-use-case/data/raw/hg38.fa.gz` to `data/spliceai/hg38.fa` and
  index it (gitignored).
- [ ] Build the 200 SNV and 120 indel cases; look up all scores; run SpliceAI on the indels.
- [ ] Results on the 200 and 120 (cross-validated). Record the model class merges.

## Phase 4: replication

- [ ] Build the 500 replication cases.
- [ ] Look up scores for the 500.
- [ ] Locked results on the 500 (primary outcome) and the checks in protocol section 9.
- [ ] EVEE descriptive analysis (protocol section 10).

## Phase 5: write-up

- [ ] `README.md` here: results, limits, and what changed from the first run and why.
- [ ] `data/README.md` and `results/README.md`: every file and field.
- [ ] Dot-and-interval charts.
- [ ] Top-level README findings, GUIDE card, and the notes in `../rare-disease/README.md` and `REPLICATION.md`;
  `docs/index.html` rebuilt.
- [ ] New release tag for the talk to cite (Tim's call).

## Phase 6 (optional, Tim's call): Evo 2 on the indel cases

- [ ] Windows for the indels only, one H100 run (about USD 3), add the scores.
- [ ] Check that no result row without Evo 2 changes (`tests/test_invariance.py`).

## Done when

- `results/` holds one result per study with the main filter, both baselines and paired intervals for every gain.
- `README.md` here describes the design, its limits, and what changed from the first run and why;
  `../rare-disease/README.md` and `REPLICATION.md` point to the results.
- The charts are dot-and-interval.
- Evo 2 can be added later: shared resamples per stratum (tested), `--fasta` required in `evo2/make_windows.py`.
