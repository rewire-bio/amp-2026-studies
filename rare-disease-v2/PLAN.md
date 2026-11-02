# Work plan

Tick each task when its output exists. Design: [`PROTOCOL.md`](PROTOCOL.md). Why: [`PROPOSAL.md`](PROPOSAL.md).

Code is copied from [`../rare-disease/`](../rare-disease/) at tag `amp-2026-v4` and changed here; each script's
docstring names the file it came from. `../rare-disease/` is not edited, so the first run stays reproducible. Large
inputs (1000 Genomes genotypes, ClinVar, MANE, the AVI table, AlphaMissense) are read from `../rare-disease/data/`,
not copied. Per-case scores already looked up there are reused; only new variants are looked up.

## Phase 0: set-up

- [x] `PROPOSAL.md`, `PROTOCOL.md`, `PLAN.md`, `README.md`
- [x] Links from the top-level `README.md`, `GUIDE.md`, `../rare-disease/README.md` and `REPLICATION.md`;
  `docs/index.html` rebuilt
- [x] `.gitignore` entries for this folder

## Phase 1: frequencies without 1000 Genomes (issue #1)

- [ ] Read the INFO header of one gnomAD v3.1.2 genomes sites VCF (`bcftools view -h` on a remote chromosome) and the
  release notes. Are there 1000 Genomes subset counts per group? Choose method A or B (protocol section 1) and log it.
- [ ] `gnomad_af.py`: AC and AN per group from the gnomAD release, minus the 1000 Genomes counts. Output
  `data/gnomad_af_no1kg.tsv` (overall AF, five group AFs, grpmax AF, found flag) for every first-run case variant
  and every ClinVar P/LP 2+ star candidate in the panel.
- [ ] Spot check: for 20 variants, compare full-gnomAD AF with the gnomAD browser, and the subtracted AF by hand.
- [ ] Redo the issue #1 table with the new frequencies: `results/gnomad_overlap.md`.

## Phase 1b: EVEE feasibility

- [ ] Check whether Zenodo serves byte ranges, so DuckDB can read the Parquet shards remotely. If not, download the
  five shards (about 35 GB) to `data/evee/` (gitignored).
- [ ] Extract EVEE rows for every first-run case variant. Report background coverage per case:
  `results/evee_coverage.md`.
- [ ] From the manuscript methods and the Zenodo record: are the released scores held-out predictions? Note it in
  `results/evee_coverage.md`.
- [ ] Apply the rule in protocol section 9; log the decision.

## Phase 2: code

- [ ] `lookup_scores.py`: sources cadd, revel, am, avi, gpnmsa, gpnstar, plus `evee` (Parquet via DuckDB, writes
  `scores_evee.tsv`); reuses rows from `../rare-disease/data/`; can score the ClinVar candidate pool.
- [ ] `build_cases.py`: panel 285 only; diagnostics chosen by CADD class (SNVs) or `annotate_indels.py` class
  (indels) from the filtered pool; background filtered by grpmax AF below 0.1% at build time; new seeds.
- [ ] `annotate_indels.py`, `run_spliceai.sh`: copied, paths updated.
- [ ] `case_benchmark.py`: one class map (protocol section 3); rule-based baseline; AVI as the primary added score
  and GPN-MSA, GPN-Star V and M as secondary; gains against both baselines; shared resamples per stratum
  (protocol section 11); `--recessive` for the 1% check; `--clinvar-covered` for the EVEE analysis; Evo 2 read
  only if score files exist.
- [ ] `avi_heldout.py`: copied, reads the new results.
- [ ] `evo2/make_windows.py`: copied with `--fasta` required.
- [ ] `figures/rd_charts.py`: dot-and-interval charts.
- [ ] Test the shared-resample property: run `case_benchmark.py` with and without one secondary score; every row
  not involving that score must be identical.

## Phase 3: build and freeze

- [ ] Look up CADD and frequencies for the whole ClinVar candidate pool (the pool CADD file is gitignored).
- [ ] Count each class in the filtered pool, with the 200's positions excluded. Do they cover the 500's quotas
  (150, 100, 75, 175)? If not, stop and log a change to the quotas before building anything.
- [ ] Build the 200 SNV and 120 indel cases: `data/cases.tsv`, `data/indel/cases.tsv`.
- [ ] Look up all scores; run SpliceAI on the indels.
- [ ] Close the **To confirm** items in protocol section 12, with Tim.
- [ ] Write `results/case_checksums.txt`. Commit and tag `rd-v2-protocol` (with Tim's approval).

## Phase 4: run

- [ ] `case_benchmark.py` on the 200 SNV cases and the 120 indel cases.
- [ ] Build and score the 500 replication cases; add their checksums.
- [ ] Locked results on the 500 (primary outcome).
- [ ] Recessive-mode check (grpmax below 1%).
- [ ] AVI held-out chromosome check.
- [ ] EVEE ClinVar-covered analysis, if Phase 1b said go.

## Phase 5: write-up

- [ ] `README.md` here: results, limits, and what changed from the first run and why.
- [ ] `data/README.md` and `results/README.md`: every file and field.
- [ ] Charts redrawn.
- [ ] Top-level `README.md` main findings and `GUIDE.md` card updated; `docs/index.html` rebuilt.
- [ ] New release tag for the talk to cite.

## Phase 6 (optional): Evo 2 on the frozen indel cases

- [ ] Apply the rule in protocol section 10.
- [ ] If yes: windows from the frozen `data/indel/case_variants.tsv`, one H100 run (about USD 11), add the scores.
- [ ] Check that every result row not involving Evo 2 is byte-for-byte unchanged.

## Open items

| Item | Where | Decided |
|---|---|---|
| Frequency method A or B | protocol section 1 | |
| Rule-based baseline scale, thresholds and tie-break | protocol section 6 | |
| EVEE coverage threshold (proposed 20) and held-out status | protocol section 9 | |
| Evo 2 threshold (proposed 80% top 1 on indels) | protocol section 10 | |

## Done when

- `results/` holds one result per study with the main filter, both baselines and paired intervals for every gain.
- `README.md` describes the design, its limits, and what changed from the first run and why.
- The charts are dot-and-interval.
- The case files are frozen and their checksums tagged; the results without Evo 2 are tagged for the talk.
