# Results

Every file in `results/`. All are committed. Charts are not here: `figures/rd_charts.py` draws them from the JSON
files into `figures/out/` (not committed), and with `--readme` copies the README's charts into `figures/` as PNG and
SVG (committed).

Every value that does not involve Evo 2 is identical to tag `amp-2026-v5` (the results before Evo 2); this version is
released as `amp-2026-v6`.

## Write-ups

| File | What it holds |
|---|---|
| [`analysis_replication.md`](analysis_replication.md) | The 500 replication cases with locked models: the primary result and verdict, the key secondary, the GPN family, Evo 2 (exploratory), groups, every section 9 check, EVEE, list lengths, ties, coverage, with p-values |
| [`analysis_snv_indel.md`](analysis_snv_indel.md) | The 200 SNV and 120 indel cases, cross-validated (descriptive), Evo 2 on both and against the first run |

## Result files (`case_benchmark.py`)

`<study>_<analysis>.json`, written by `case_benchmark.py --study <study> --analysis <analysis>`.

| Study | Analyses | Models |
|---|---|---|
| `snv` (200 cases, first run's 50 people) | main, bg1pct, nomoi, all1kg, heldout, evee | 5-fold cross-validation grouped by person |
| `indel` (120 cases, same people) | main, bg1pct, nomoi, all1kg | 5-fold cross-validation grouped by person; SpliceAI, Evo 2, Rule-based, Lookups combined and Lookups + Evo 2 |
| `replication` (500 cases, 50 new people) | main, bg1pct, nomoi, all1kg, heldout, evee | Locked (fitted on the 200 SNV cases) and, as `Lookups combined` without "(locked)", cross-validated on the 500 |

| Analysis | Filters |
|---|---|
| `main` | Background and diagnostic grpmax below 0.1% (background with the person and family left out of gnomAD); background inheritance filter |
| `bg1pct` | grpmax below 1%; locked models fitted under 1% |
| `nomoi` | grpmax below 0.1%; no inheritance filter |
| `all1kg` | Background grpmax with all of 1000 Genomes left out below 0.1%; inheritance filter |
| `heldout` | Main filters; each case ranked only among background variants on chromosomes of its diagnostic's split (AVI training chromosomes 1, 4, 7, 8, 10, 13, 15, or the others) |
| `evee` | bg1pct filters; lists restricted to variants with an EVEE row; cases whose diagnostic has none dropped. Descriptive only |

Fields in each JSON. Proportions are fractions (0.906 is 90.6%); intervals are `[low, high]`.

| Field | Contents |
|---|---|
| `study`, `analysis`, `cases`, `people`, `median_variants_per_case` | What was run; median list length, diagnostic included |
| `filter` | `max_af`, the frequency columns used, whether the inheritance filter ran, and row counts at each step: `variants_before`, `variants_after_frequency`, `variants_after_inheritance`, the other type's rows counted for inheritance (`background_indels_in_inheritance_count` or `background_snvs_in_inheritance_count`), `variants_after_dropping_other_type`; `cases_lost_diagnostic` (0 everywhere) |
| `load_notes` | `case_variants_without_mane_cadd_row`, `diagnostic_class_differs_from_cases_tsv` (0 everywhere) |
| `coverage` | Per score: `variants` and `diagnostics` with a score in the ranked lists |
| `scores_without_coverage`, `dropped_scores` | Scores with no rows for the study (the indel study: all but SpliceAI and Evo 2); scores removed with `--drop-score` |
| `model_classes` | Diagnostics per class in the training set, classes merged into other, the model classes |
| `fits`, `intervals`, `package_versions` | How models were fitted; the bootstrap seed rule; package versions |
| `rule_ties` | Not in heldout. Share of cases with a background variant TIED with the diagnostic, for Rule-based and Rule + AVI (empty in evee, which has no rule scores) |
| `strata` | Not in heldout. One entry per stratum: `all`, each diagnostic group, each CADD class within non-coding (`noncoding: <class>`), and `gnomad: found` and `gnomad: absent`. Each has `cases`, `people`, `small` (under 20 cases), `scores`, `rule_ties`, `seed_label` and, except in the gnomAD strata and in evee, `contrasts`. evee has no gnomAD strata |
| `strata.<s>.scores.<score>` | `top_1`, `top_5`, `top_10` with `_ci95`, and `median_competitors` |
| `strata.<s>.contrasts.<A minus B>` | Paired gains `top_1`, `top_5`, `top_10` with `_ci95` and `_p` (two-sided bootstrap p-value on the same draws); `role` (primary, key secondary, secondary, descriptive, or for Evo 2 "exploratory, Bonferroni-adjusted over four scores"); the GPN family also has `_ci98.33`, Evo 2 `_ci98.75`; the locked primary has `confirmatory` and `verdict` |
| `primary` | Replication only: the locked primary contrast with its `verdict` (the reading rule of PROTOCOL.md section 8; it applies to `main`) |
| `gnomad_split` | Not in heldout or evee. `diagnostics_found` (in full gnomAD with AC above 0), `diagnostics_absent`; the primary and key secondary in each stratum on one shared draw, with `_p`; `difference` (absent minus found) with `top_1_p` and `possibly_gnomad_presence` |
| `heldout` | heldout only: `training_chromosomes` and `other_chromosomes`, each with `cases`, `median_variants_per_case`, `scores`, `contrasts`; `difference` (training minus other) with `top_1_p` and `possibly_inflated_by_training_overlap` |
| `evee_note`, `cases_without_evee_diagnostic` | evee only |

Score names: CADD, SpliceAI, REVEL, AlphaMissense (the lookups); AVI, GPN-MSA, GPN-Star-V, GPN-Star-M, Evo 2 (added
scores; Evo 2 is minus `evo2_delta` from `../evo2/run*/scores*.tsv`, higher is more damaging); EVEE (evee only);
Rule-based and Rule + AVI; Lookups combined and Lookups + X; on the replication study also the same models with
"(locked)".

Levels: the primary is read at 95% by the rule in PROTOCOL.md section 8; the GPN family at 98.33%; Evo 2 at 98.75%,
exploratory (section 12); everything else at 95%, descriptive. The 16 files hold 944 paired top-1 contrasts
(`strata`, `gnomad_split` and `heldout` together), 175 of them with p below 0.05.

## Build and check reports

| File | Made by | What it holds |
|---|---|---|
| `build_snv.md`, `build_indel.md`, `build_replication.md` | `build_cases.py` | Cases, people, seed, background per case before and after the main filter, quotas and draws per group, inheritance and zygosity counts, positions excluded (replication), and the build checks |
| `pool_counts.md` | `pool.py` | Eligible ClinVar variants and positions per group against the quotas, the non-coding group by CADD class, groups by inheritance mode, and ClinVar's consequence against this run's class |
| `gnomad_checks.md` | `gnomad_af.py checks` | The section 1 checks: adj counts against the subset file's own AC and AN, full counts against MyVariant.info, released 1000 Genomes counts within gnomAD's group counts, group sizes for the all-1000-Genomes fallback, matching and clamping, and hand-worked examples |
| `gnomad_overlap.md` | `gnomad_af.py overlap` | Issue #1's table redone with the new frequencies, and background per case after the grpmax filter by frequency source |
| `evee_source.md` | `lookup_scores.py evee` | The EVEE Zenodo record, download date, row counts and checksums |

## Logs

Run logs are in `../logs/`: `bench_<study>_<analysis>.log` for each result file, `cb_invariance.log` for
`tests/test_invariance.py` (run before Evo 2 was added), `evo2_launch.log` for the Evo 2 run, and one log per build
and lookup step.
