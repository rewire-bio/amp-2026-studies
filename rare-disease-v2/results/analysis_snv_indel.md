# Results on the 200 SNV and 120 indel cases (cross-validated)

These are the 200 SNV and 120 indel cases built from the first run's 50 people, with every model fitted by 5-fold cross-validation grouped by person. They are descriptive. The primary outcome is the locked models on the 500 replication cases, in [`analysis_replication.md`](analysis_replication.md) (PROTOCOL.md section 8). No verdict is read off these numbers.

Made by `case_benchmark.py` on the committed case and score files: `results/snv_<analysis>.json` for main, bg1pct, nomoi, all1kg, heldout and evee, and `results/indel_<analysis>.json` for main, bg1pct, nomoi and all1kg. Logs: `logs/bench_<study>_<analysis>.log`.

All figures are percentages, rounded to one decimal; gains are in percentage points. Intervals are 95% bootstrap intervals over people (1,000 resamples, shared by every score in a stratum); every gain has a paired interval and a two-sided bootstrap p-value (`top_k_p`). The GPN family also has a Bonferroni-adjusted 98.33% interval, and Evo 2, exploratory (PROTOCOL.md section 12), a 98.75% interval. Everything here is descriptive: `snv_main.json` holds 78 paired top-1 contrasts and `indel_main.json` 14, so some intervals will exclude zero by chance. Every value not involving Evo 2 is identical to tag `amp-2026-v5`. Top k is the share of cases with fewer than k competitors. For single scores and models a competitor is a background variant scoring at least as high as the diagnostic; for the rule scores it is a background variant whose pairwise comparison with the diagnostic is TIED or BEATS (section 7). Top 1 is the share of cases where the diagnostic ranks first.

In the indel study, "SpliceAI" is the score we ran, and each indel diagnostic is ranked among the person's background indels only. CADD and precomputed SpliceAI are set aside for indels by design (section 6), and REVEL, AlphaMissense, AVI and the GPN scores have no indel rows. So the indel study has SpliceAI, Evo 2, Rule-based, Lookups combined and Lookups + Evo 2, one added-score contrast (Evo 2), and no AVI contrast.

Cross-validated numbers move by a few points when people move between folds. `GroupKFold` assigns people to folds by row count, so a small change to the lists (a rebuild that adds a few background rows) can move many people to another fold and change Lookups combined by a few cases. This is one reason these numbers are descriptive.

## Headline (main analysis)

- SNV, 200 cases, median list 26 variants (diagnostic included). Top 1: Lookups combined 84.0 (75.5 to 91.5), Rule-based 81.0 (73.5 to 88.0), AVI alone 79.0 (71.0 to 86.0), Lookups + AVI 88.0 (81.0 to 94.0), Rule + AVI 85.5 (79.0 to 91.0).
- Primary contrast, cross-validated: Lookups + AVI minus Lookups combined, top 1 +4.0 (-2.0 to +9.5), p = 0.178.
- Key secondary: Rule + AVI minus Rule-based, top 1 +4.5 (+0.5 to +8.0), p = 0.028.
- GPN family added to Lookups combined (98.33% interval): GPN-MSA -2.0 (-8.0 to +3.5), GPN-Star-V -0.5 (-6.0 to +5.0), GPN-Star-M +0.5 (-4.5 to +5.5). None excludes zero.
- Evo 2 added to Lookups combined (exploratory, 98.75% interval): SNV +1.0 (-4.5 to +6.0), p = 0.760; indel +1.7 (-9.0 to +12.8), p = 0.746. Evo 2 alone: SNV top 1 51.5 (43.5 to 59.5), indel 81.7 (73.7 to 89.2).
- The key secondary gain sits in the cases whose diagnostic is absent from gnomAD (+6.1 (+1.9 to +11.5); 0.0 (-6.7 to +5.0) in found cases). See the gnomAD split below.
- Indel, 120 cases, median list 4 indels. Top 1: Lookups combined 84.2 (75.4 to 92.4), Rule-based 76.7 (65.8 to 86.8), SpliceAI 58.3 (49.6 to 66.4). Lookups combined minus Rule-based +7.5 (+3.2 to +12.6), p = 0.001.

## Cases and list lengths

Median list length per case, diagnostic included. No case lost its diagnostic under any filter (`cases_lost_diagnostic` is 0 in every run).

| Analysis | SNV cases | SNV median list | SNV variants in lists | Indel cases | Indel median list | Indel variants in lists |
|---|---|---|---|---|---|---|
| main | 200 | 26 | 5,401 | 120 | 4 | 576 |
| bg1pct | 200 | 78.5 | 19,561 | 120 | 14 | 1,955 |
| nomoi | 200 | 51.5 | 9,816 | 120 | 7 | 781 |
| all1kg | 200 | 25.5 | 5,581 | 120 | 4 | 610 |
| heldout | 200 | 7 (training split), 19 (other) | 5,401 | not run (SNV only) | | |
| evee | 200 | 29.5 | 7,084 (with an EVEE row) | not run (SNV only) | | |

The inheritance filter counts both variant types, then keeps only the study's own type. SNV main: 10,956 rows after the frequency filter, 6,189 after the inheritance filter (788 of them background indels), 5,401 SNVs ranked. Indel main: 6,410, then 3,615 (3,039 background SNVs), 576 indels ranked.

## Model class merges

- SNV (from the 200): canonical splice 40, intronic 38, missense 60, nonsense 30, splice region 23, synonymous 6, utr 3. Merged into other: synonymous, utr. Model classes: canonical splice, intronic, missense, nonsense, other, splice region.
- Indel (from the 120): canonical splice 30, frameshift 40, inframe 30, intronic 13, utr 7. Merged into other: utr. Model classes: canonical splice, frameshift, inframe, intronic, other.
- The same merges hold in every analysis. `diagnostic_class_differs_from_cases_tsv` and `case_variants_without_mane_cadd_row` are 0 in every run.

## Rule tie share

Share of cases with at least one background variant whose pairwise outcome is TIED.

| Analysis | SNV Rule-based | SNV Rule + AVI | Indel Rule-based |
|---|---|---|---|
| main | 0.0% (0 of 200) | 0.0% (0 of 200) | 4.2% (5 of 120) |
| bg1pct | 0.0% (0 of 200) | 0.0% (0 of 200) | 9.2% (11 of 120) |
| nomoi | 1.0% (2 of 200) | 0.0% (0 of 200) | 5.0% (6 of 120) |
| all1kg | 0.0% (0 of 200) | 0.0% (0 of 200) | 4.2% (5 of 120) |

Indel ties in main by group: frameshift 5.0% (2 of 40), inframe 10.0% (3 of 30), canonical splice 0.0% (0 of 30), non-coding 0.0% (0 of 20). They are indels with the same points and the same SpliceAI score (often 0) as the diagnostic; indels have no CADD to break them.

## Coverage per score (main analysis)

Variants with a score over all variants in the lists (diagnostics included), and diagnostics with a score. CADD is not used for indels.

| Score | SNV variants (of 5,401) | SNV diagnostics (of 200) | Indel variants (of 576) | Indel diagnostics (of 120) |
|---|---|---|---|---|
| CADD | 5,401 | 200 | not used | not used |
| SpliceAI | 4,485 | 192 | 483 | 111 |
| REVEL | 1,244 | 68 | none | none |
| AlphaMissense | 1,025 | 57 | none | none |
| AVI | 5,401 | 200 | none | none |
| GPN-MSA | 5,401 | 200 | none | none |
| GPN-Star-V | 5,401 | 200 | none | none |
| GPN-Star-M | 5,401 | 200 | none | none |
| Evo 2 | 5,401 | 200 | 576 | 120 |
| EVEE (evee analysis) | 7,084 (the EVEE-covered bg1pct lists, of 19,561) | 200 | not used | not used |

9 of the 120 indel diagnostics have no SpliceAI score, so SpliceAI and the rule-based ranking put them below every scored indel. Evo 2 scores every SNV and indel in every analysis.

## Results by analysis, all cases

### SNV study, 200 cases

#### main (grpmax below 0.1%, family left out, inheritance filter)

200 cases, median list 26.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 51.0 (43.0 to 59.5) | 82.0 (76.0 to 87.5) | 93.0 (90.0 to 96.0) |
| SpliceAI | 53.5 (46.0 to 59.5) | 67.0 (61.0 to 72.0) | 70.5 (64.0 to 76.5) |
| REVEL | 27.5 (23.0 to 32.5) | 32.5 (27.5 to 37.5) | 38.0 (31.5 to 45.5) |
| AlphaMissense | 23.5 (19.0 to 28.0) | 28.5 (24.5 to 33.5) | 32.5 (26.5 to 39.5) |
| AVI | 79.0 (71.0 to 86.0) | 95.0 (92.0 to 98.0) | 99.0 (97.5 to 100.0) |
| GPN-MSA | 45.5 (36.0 to 55.0) | 76.5 (70.5 to 82.5) | 88.5 (83.5 to 93.0) |
| GPN-Star-V | 52.5 (44.0 to 60.5) | 85.5 (80.0 to 90.0) | 92.0 (88.0 to 96.0) |
| GPN-Star-M | 50.5 (42.5 to 58.0) | 80.5 (74.0 to 86.5) | 93.5 (89.5 to 96.5) |
| Evo 2 | 51.5 (43.5 to 59.5) | 82.0 (76.0 to 88.0) | 94.5 (90.5 to 97.5) |
| Rule-based | 81.0 (73.5 to 88.0) | 96.0 (92.5 to 99.0) | 97.5 (94.5 to 99.5) |
| Rule + AVI | 85.5 (79.0 to 91.0) | 96.5 (93.5 to 99.0) | 97.5 (94.5 to 99.5) |
| Lookups combined | 84.0 (75.5 to 91.5) | 97.0 (94.5 to 99.5) | 98.5 (96.0 to 100.0) |
| Lookups + AVI | 88.0 (81.0 to 94.0) | 98.0 (95.5 to 100.0) | 98.0 (95.5 to 100.0) |
| Lookups + GPN-MSA | 82.0 (74.0 to 89.0) | 96.0 (93.0 to 98.5) | 97.0 (94.0 to 99.5) |
| Lookups + GPN-Star-V | 83.5 (75.5 to 91.0) | 97.0 (94.5 to 99.5) | 97.5 (95.0 to 99.5) |
| Lookups + GPN-Star-M | 84.5 (75.5 to 91.5) | 97.0 (94.0 to 99.5) | 98.0 (95.5 to 100.0) |
| Lookups + Evo 2 | 85.0 (77.0 to 91.5) | 97.0 (94.0 to 99.5) | 98.5 (96.0 to 100.0) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Top 1, adjusted interval |
|---|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast, cross-validated (descriptive) | +4.0 (-2.0 to +9.5), p = 0.178 | +1.0 (-1.0 to +2.5), p = 0.454 | -0.5 (-1.5 to 0.0), p = 0.726 |  |
| Rule + AVI minus Rule-based | key secondary | +4.5 (+0.5 to +8.0), p = 0.028 | +0.5 (-1.0 to +2.5), p = 0.754 | 0.0 (0.0 to 0.0), p = 1.000 |  |
| Lookups combined minus Rule-based | descriptive | +3.0 (-4.5 to +9.5), p = 0.434 | +1.0 (-1.5 to +4.0), p = 0.568 | +1.0 (0.0 to +2.5), p = 0.290 |  |
| Lookups + AVI minus Rule-based | descriptive | +7.0 (0.0 to +13.0), p = 0.052 | +2.0 (-0.5 to +5.0), p = 0.206 | +0.5 (-1.0 to +2.0), p = 0.802 |  |
| Lookups + GPN-MSA minus Lookups combined | secondary, cross-validated (descriptive) | -2.0 (-7.0 to +2.5), p = 0.476 | -1.0 (-2.5 to 0.0), p = 0.294 | -1.5 (-3.5 to 0.0), p = 0.094 | 98.33% (-8.0 to +3.5) |
| Lookups + GPN-Star-V minus Lookups combined | secondary, cross-validated (descriptive) | -0.5 (-5.5 to +4.0), p = 0.934 | 0.0 (-1.5 to +1.5), p = 1.000 | -1.0 (-2.5 to 0.0), p = 0.290 | 98.33% (-6.0 to +5.0) |
| Lookups + GPN-Star-M minus Lookups combined | secondary, cross-validated (descriptive) | +0.5 (-4.0 to +4.5), p = 0.958 | 0.0 (-1.5 to +1.5), p = 1.000 | -0.5 (-1.5 to 0.0), p = 0.742 | 98.33% (-4.5 to +5.5) |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | +1.0 (-3.5 to +5.5), p = 0.760 | 0.0 (-1.5 to +1.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 | 98.75% (-4.5 to +6.0) |

#### bg1pct (grpmax below 1%, inheritance filter)

200 cases, median list 78.5.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 41.0 (33.5 to 48.5) | 67.0 (60.0 to 74.0) | 80.0 (73.0 to 86.5) |
| SpliceAI | 47.0 (40.0 to 54.5) | 61.5 (56.0 to 67.0) | 64.0 (58.5 to 69.5) |
| REVEL | 25.5 (20.5 to 30.5) | 30.5 (25.5 to 36.0) | 32.0 (27.0 to 37.5) |
| AlphaMissense | 22.5 (18.0 to 27.5) | 28.0 (23.5 to 33.0) | 28.5 (24.0 to 33.5) |
| AVI | 70.0 (61.5 to 78.0) | 93.0 (89.0 to 96.5) | 95.5 (92.0 to 98.5) |
| GPN-MSA | 34.5 (27.5 to 42.5) | 62.0 (55.5 to 68.0) | 76.0 (72.0 to 80.0) |
| GPN-Star-V | 45.0 (38.0 to 52.0) | 71.0 (65.0 to 77.0) | 84.5 (80.0 to 89.0) |
| GPN-Star-M | 38.0 (30.0 to 46.0) | 68.0 (62.0 to 73.5) | 83.5 (79.0 to 88.0) |
| Evo 2 | 39.0 (32.5 to 46.0) | 65.5 (59.0 to 72.0) | 80.5 (74.0 to 86.5) |
| Rule-based | 71.5 (62.5 to 79.5) | 93.5 (89.5 to 97.0) | 96.5 (93.5 to 99.0) |
| Rule + AVI | 75.5 (68.0 to 82.5) | 95.5 (92.0 to 98.0) | 96.0 (93.0 to 98.5) |
| Lookups combined | 77.0 (68.0 to 85.0) | 95.5 (92.0 to 98.5) | 96.5 (93.5 to 99.0) |
| Lookups + AVI | 78.0 (69.5 to 85.5) | 95.5 (92.0 to 98.0) | 96.0 (92.5 to 98.5) |
| Lookups + GPN-MSA | 76.5 (68.5 to 84.0) | 95.5 (92.0 to 98.0) | 96.0 (92.5 to 98.5) |
| Lookups + GPN-Star-V | 77.0 (69.5 to 84.0) | 95.0 (91.0 to 98.0) | 95.0 (91.0 to 98.0) |
| Lookups + GPN-Star-M | 77.5 (70.5 to 84.0) | 94.5 (90.5 to 97.5) | 95.5 (92.0 to 98.0) |
| Lookups + Evo 2 | 77.5 (70.0 to 85.0) | 95.0 (91.0 to 98.0) | 95.5 (92.0 to 98.5) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Top 1, adjusted interval |
|---|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast, cross-validated (descriptive) | +1.0 (-7.5 to +9.5), p = 0.850 | 0.0 (-2.0 to +2.0), p = 1.000 | -0.5 (-3.0 to +1.5), p = 0.844 |  |
| Rule + AVI minus Rule-based | key secondary | +4.0 (+0.5 to +8.0), p = 0.038 | +2.0 (0.0 to +4.5), p = 0.124 | -0.5 (-1.5 to 0.0), p = 0.734 |  |
| Lookups combined minus Rule-based | descriptive | +5.5 (-2.5 to +13.5), p = 0.212 | +2.0 (-1.0 to +5.5), p = 0.300 | 0.0 (-2.0 to +2.0), p = 1.000 |  |
| Lookups + AVI minus Rule-based | descriptive | +6.5 (-2.5 to +15.5), p = 0.168 | +2.0 (-1.0 to +5.5), p = 0.240 | -0.5 (-2.0 to +1.0), p = 0.774 |  |
| Lookups + GPN-MSA minus Lookups combined | secondary, cross-validated (descriptive) | -0.5 (-5.0 to +4.0), p = 0.920 | 0.0 (-1.5 to +1.5), p = 1.000 | -0.5 (-1.5 to 0.0), p = 0.734 | 98.33% (-6.5 to +5.5) |
| Lookups + GPN-Star-V minus Lookups combined | secondary, cross-validated (descriptive) | 0.0 (-6.0 to +6.0), p = 1.000 | -0.5 (-2.5 to +1.0), p = 0.774 | -1.5 (-3.5 to 0.0), p = 0.108 | 98.33% (-7.0 to +7.0) |
| Lookups + GPN-Star-M minus Lookups combined | secondary, cross-validated (descriptive) | +0.5 (-4.5 to +6.0), p = 0.978 | -1.0 (-2.5 to 0.0), p = 0.254 | -1.0 (-3.0 to +1.0), p = 0.410 | 98.33% (-5.5 to +7.5) |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | +0.5 (-3.5 to +5.5), p = 0.970 | -0.5 (-1.5 to 0.0), p = 0.722 | -1.0 (-3.0 to +0.5), p = 0.442 | 98.75% (-4.5 to +7.0) |

#### nomoi (grpmax below 0.1%, no inheritance filter)

200 cases, median list 51.5.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 40.0 (32.5 to 47.0) | 73.5 (66.5 to 80.0) | 83.5 (78.0 to 89.0) |
| SpliceAI | 46.0 (39.0 to 53.0) | 63.0 (57.5 to 68.5) | 66.5 (61.0 to 72.0) |
| REVEL | 23.5 (18.0 to 29.0) | 30.5 (25.5 to 36.0) | 33.0 (28.0 to 38.5) |
| AlphaMissense | 20.0 (15.0 to 25.0) | 28.0 (23.5 to 32.5) | 28.0 (23.5 to 32.5) |
| AVI | 66.5 (58.5 to 74.0) | 91.0 (86.0 to 95.0) | 95.0 (92.0 to 98.0) |
| GPN-MSA | 35.5 (27.0 to 44.0) | 67.5 (61.0 to 74.0) | 82.0 (76.5 to 87.5) |
| GPN-Star-V | 44.0 (36.0 to 51.5) | 77.5 (71.5 to 83.5) | 87.0 (82.5 to 91.0) |
| GPN-Star-M | 39.5 (32.0 to 47.0) | 72.0 (65.0 to 78.5) | 83.0 (77.5 to 88.5) |
| Evo 2 | 39.5 (31.0 to 48.5) | 72.5 (65.5 to 79.0) | 83.5 (78.0 to 88.5) |
| Rule-based | 66.5 (58.0 to 75.0) | 93.0 (89.0 to 96.5) | 96.5 (93.0 to 99.0) |
| Rule + AVI | 73.5 (65.5 to 81.5) | 94.0 (90.5 to 97.5) | 96.0 (93.0 to 98.5) |
| Lookups combined | 68.5 (57.5 to 79.0) | 95.0 (91.5 to 98.0) | 96.5 (93.5 to 99.0) |
| Lookups + AVI | 74.0 (65.0 to 82.5) | 94.5 (91.0 to 97.5) | 96.5 (93.5 to 99.0) |
| Lookups + GPN-MSA | 71.0 (61.5 to 80.5) | 93.5 (89.5 to 97.5) | 95.5 (92.5 to 98.5) |
| Lookups + GPN-Star-V | 71.0 (62.0 to 79.5) | 94.5 (91.0 to 97.5) | 96.5 (94.0 to 98.5) |
| Lookups + GPN-Star-M | 68.0 (59.0 to 76.5) | 93.0 (89.0 to 96.5) | 95.5 (92.5 to 98.5) |
| Lookups + Evo 2 | 74.0 (65.0 to 83.0) | 93.0 (89.0 to 96.5) | 95.0 (91.5 to 98.0) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Top 1, adjusted interval |
|---|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast, cross-validated (descriptive) | +5.5 (-2.5 to +13.5), p = 0.210 | -0.5 (-2.0 to +1.0), p = 0.776 | 0.0 (-2.0 to +2.0), p = 1.000 |  |
| Rule + AVI minus Rule-based | key secondary | +7.0 (+2.5 to +12.0), p = 0.002 | +1.0 (0.0 to +2.5), p = 0.282 | -0.5 (-2.0 to +1.0), p = 0.804 |  |
| Lookups combined minus Rule-based | descriptive | +2.0 (-7.0 to +12.5), p = 0.708 | +2.0 (-0.5 to +5.0), p = 0.186 | 0.0 (-2.0 to +2.5), p = 1.000 |  |
| Lookups + AVI minus Rule-based | descriptive | +7.5 (-1.5 to +16.5), p = 0.110 | +1.5 (-1.5 to +5.0), p = 0.428 | 0.0 (-2.5 to +3.0), p = 1.000 |  |
| Lookups + GPN-MSA minus Lookups combined | secondary, cross-validated (descriptive) | +2.5 (-1.5 to +7.0), p = 0.280 | -1.5 (-3.5 to 0.0), p = 0.082 | -1.0 (-2.5 to 0.0), p = 0.252 | 98.33% (-2.0 to +8.0) |
| Lookups + GPN-Star-V minus Lookups combined | secondary, cross-validated (descriptive) | +2.5 (-3.0 to +8.0), p = 0.408 | -0.5 (-1.5 to 0.0), p = 0.700 | 0.0 (-1.5 to +1.5), p = 1.000 | 98.33% (-4.0 to +9.0) |
| Lookups + GPN-Star-M minus Lookups combined | secondary, cross-validated (descriptive) | -0.5 (-5.5 to +4.5), p = 0.934 | -2.0 (-4.5 to 0.0), p = 0.122 | -1.0 (-2.5 to 0.0), p = 0.252 | 98.33% (-6.8 to +5.5) |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | +5.5 (+1.0 to +10.5), p = 0.014 | -2.0 (-4.0 to -0.5), p = 0.040 | -1.5 (-3.5 to 0.0), p = 0.118 | 98.75% (+0.1 to +12.0) |

#### all1kg (all of 1000 Genomes left out of the background frequency, inheritance filter)

200 cases, median list 25.5.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 48.5 (40.5 to 56.5) | 81.5 (75.0 to 87.5) | 93.0 (89.5 to 96.0) |
| SpliceAI | 52.5 (45.5 to 58.0) | 67.0 (61.0 to 72.0) | 68.0 (62.0 to 73.5) |
| REVEL | 27.5 (23.0 to 33.0) | 32.5 (27.5 to 38.0) | 36.0 (30.0 to 43.5) |
| AlphaMissense | 23.5 (19.5 to 28.5) | 28.5 (24.5 to 33.0) | 31.0 (25.0 to 38.0) |
| AVI | 77.5 (69.0 to 85.5) | 95.0 (91.5 to 98.0) | 99.0 (97.5 to 100.0) |
| GPN-MSA | 43.5 (35.0 to 52.5) | 76.5 (71.0 to 82.5) | 88.5 (83.5 to 93.0) |
| GPN-Star-V | 52.5 (44.5 to 60.5) | 85.0 (80.0 to 89.5) | 91.5 (87.0 to 95.5) |
| GPN-Star-M | 50.0 (42.5 to 57.0) | 80.0 (73.5 to 85.5) | 93.0 (89.0 to 96.5) |
| Evo 2 | 52.0 (44.0 to 60.0) | 82.0 (76.5 to 87.5) | 94.0 (90.5 to 97.5) |
| Rule-based | 79.5 (70.5 to 86.5) | 96.0 (92.5 to 99.0) | 97.5 (94.5 to 99.5) |
| Rule + AVI | 84.0 (77.0 to 90.0) | 96.5 (93.5 to 99.0) | 97.5 (94.5 to 99.5) |
| Lookups combined | 77.0 (68.0 to 85.5) | 97.5 (95.0 to 99.5) | 98.5 (96.0 to 100.0) |
| Lookups + AVI | 84.5 (76.5 to 91.0) | 97.5 (94.5 to 99.5) | 97.5 (94.5 to 99.5) |
| Lookups + GPN-MSA | 79.0 (70.5 to 86.5) | 97.0 (94.0 to 99.5) | 98.0 (95.5 to 100.0) |
| Lookups + GPN-Star-V | 82.5 (74.5 to 89.5) | 98.0 (95.5 to 100.0) | 99.0 (97.5 to 100.0) |
| Lookups + GPN-Star-M | 82.5 (75.0 to 89.5) | 97.0 (94.0 to 99.5) | 98.0 (96.0 to 99.5) |
| Lookups + Evo 2 | 76.0 (66.5 to 85.0) | 97.0 (94.5 to 99.5) | 98.5 (96.0 to 100.0) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Top 1, adjusted interval |
|---|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast, cross-validated (descriptive) | +7.5 (+2.0 to +14.5), p = 0.012 | 0.0 (-2.0 to +2.0), p = 1.000 | -1.0 (-2.5 to 0.0), p = 0.270 |  |
| Rule + AVI minus Rule-based | key secondary | +4.5 (+1.0 to +8.5), p = 0.022 | +0.5 (-1.0 to +2.0), p = 0.736 | 0.0 (0.0 to 0.0), p = 1.000 |  |
| Lookups combined minus Rule-based | descriptive | -2.5 (-11.0 to +5.5), p = 0.602 | +1.5 (-1.0 to +4.5), p = 0.322 | +1.0 (0.0 to +2.5), p = 0.238 |  |
| Lookups + AVI minus Rule-based | descriptive | +5.0 (-2.5 to +12.0), p = 0.196 | +1.5 (-1.0 to +4.5), p = 0.356 | 0.0 (-2.0 to +2.0), p = 1.000 |  |
| Lookups + GPN-MSA minus Lookups combined | secondary, cross-validated (descriptive) | +2.0 (-2.0 to +6.0), p = 0.422 | -0.5 (-1.5 to 0.0), p = 0.668 | -0.5 (-1.5 to 0.0), p = 0.668 | 98.33% (-2.5 to +7.0) |
| Lookups + GPN-Star-V minus Lookups combined | secondary, cross-validated (descriptive) | +5.5 (+0.5 to +11.0), p = 0.038 | +0.5 (0.0 to +1.5), p = 0.748 | +0.5 (0.0 to +1.5), p = 0.762 | 98.33% (-0.3 to +12.5) |
| Lookups + GPN-Star-M minus Lookups combined | secondary, cross-validated (descriptive) | +5.5 (+1.5 to +9.5), p = 0.001 | -0.5 (-1.5 to 0.0), p = 0.668 | -0.5 (-2.0 to +1.0), p = 0.676 | 98.33% (+1.0 to +10.5) |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | -1.0 (-5.5 to +3.5), p = 0.758 | -0.5 (-1.5 to 0.0), p = 0.804 | 0.0 (0.0 to 0.0), p = 1.000 | 98.75% (-7.0 to +4.5) |

The adjusted column is the Bonferroni 98.33% interval for the GPN family (section 8) and the 98.75% interval for Evo 2 (section 12). Two exclude zero. Evo 2 in nomoi: +5.5 (+1.0 to +10.5), p = 0.014, 98.75% (+0.1 to +12.0); top 5 moves the other way (-2.0, -4.0 to -0.5), and on the 500 under the same check Evo 2 adds 0.0 (-4.6 to +3.2), 98.75% (-6.1 to +3.8). The other is GPN-Star-M in all1kg, +5.5 (+1.5 to +9.5), p = 0.001, 98.33% (+1.0 to +10.5). It is a cross-validated check on the 200; on the 500 under the same check GPN-Star-M adds nothing (`analysis_replication.md`).

### Indel study, 120 cases, ranked among background indels

Rule-based top 5 can be below 100% although the median list is short: a diagnostic whose class needs SpliceAI and has none ranks below every scored indel.

#### main (grpmax below 0.1%, family left out, inheritance filter)

120 cases, median list 4.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 58.3 (49.6 to 66.4) | 90.8 (84.5 to 96.0) | 99.2 (97.4 to 100.0) |
| Evo 2 | 81.7 (73.7 to 89.2) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) |
| Rule-based | 76.7 (65.8 to 86.8) | 95.8 (91.8 to 99.2) | 99.2 (97.4 to 100.0) |
| Lookups combined | 84.2 (75.4 to 92.4) | 95.8 (91.8 to 99.2) | 99.2 (97.4 to 100.0) |
| Lookups + Evo 2 | 85.8 (77.7 to 93.5) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Adjusted interval |
|---|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +7.5 (+3.2 to +12.6), p = 0.001 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 |  |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | +1.7 (-5.9 to +9.9), p = 0.746 | +4.2 (+0.8 to +8.2), p = 0.012 | +0.8 (0.0 to +2.6), p = 0.768 | 98.75% (-9.0 to +12.8); top 5 (+0.8 to +9.5) |

#### bg1pct (grpmax below 1%, inheritance filter)

120 cases, median list 14.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 42.5 (33.6 to 51.7) | 65.8 (57.7 to 74.6) | 75.8 (66.7 to 83.6) |
| Evo 2 | 75.0 (67.5 to 81.7) | 98.3 (95.9 to 100.0) | 99.2 (97.4 to 100.0) |
| Rule-based | 70.8 (59.8 to 80.0) | 94.2 (89.8 to 97.6) | 94.2 (89.8 to 97.6) |
| Lookups combined | 75.0 (65.5 to 83.6) | 94.2 (89.8 to 97.6) | 94.2 (89.8 to 97.6) |
| Lookups + Evo 2 | 82.5 (73.8 to 90.0) | 99.2 (97.4 to 100.0) | 99.2 (97.4 to 100.0) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Adjusted interval |
|---|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +4.2 (-1.7 to +10.4), p = 0.236 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 |  |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | +7.5 (0.0 to +15.8), p = 0.064 | +5.0 (+0.8 to +9.8), p = 0.024 | +5.0 (+0.8 to +9.8), p = 0.024 | 98.75% (-2.4 to +18.1); top 5 (0.0 to +11.6) |

#### nomoi (grpmax below 0.1%, no inheritance filter)

120 cases, median list 7.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 53.3 (44.9 to 62.0) | 79.2 (69.9 to 87.4) | 99.2 (97.4 to 100.0) |
| Evo 2 | 78.3 (70.8 to 86.2) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) |
| Rule-based | 74.2 (63.6 to 84.5) | 95.0 (91.0 to 98.4) | 99.2 (97.4 to 100.0) |
| Lookups combined | 81.7 (73.0 to 89.9) | 95.0 (91.0 to 98.4) | 99.2 (97.4 to 100.0) |
| Lookups + Evo 2 | 87.5 (81.2 to 93.3) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Adjusted interval |
|---|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +7.5 (+2.5 to +12.9), p = 0.004 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 |  |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | +5.8 (0.0 to +12.1), p = 0.058 | +5.0 (+1.6 to +9.0), p = 0.008 | +0.8 (0.0 to +2.6), p = 0.738 | 98.75% (-0.9 to +13.8); top 5 (+0.8 to +10.2) |

#### all1kg (all of 1000 Genomes left out of the background frequency, inheritance filter)

120 cases, median list 4.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 59.2 (50.8 to 68.6) | 89.2 (82.3 to 95.7) | 99.2 (97.4 to 100.0) |
| Evo 2 | 81.7 (73.6 to 89.1) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) |
| Rule-based | 76.7 (66.1 to 86.9) | 95.8 (92.2 to 99.2) | 99.2 (97.4 to 100.0) |
| Lookups combined | 83.3 (75.4 to 90.8) | 95.8 (92.2 to 99.2) | 99.2 (97.4 to 100.0) |
| Lookups + Evo 2 | 88.3 (80.8 to 94.7) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) |

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain | Adjusted interval |
|---|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +6.7 (+2.5 to +11.5), p = 0.002 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 |  |
| Lookups + Evo 2 minus Lookups combined | exploratory, cross-validated | +5.0 (-2.5 to +12.1), p = 0.224 | +4.2 (+0.8 to +7.8), p = 0.010 | +0.8 (0.0 to +2.6), p = 0.750 | 98.75% (-4.1 to +14.2); top 5 (+0.8 to +9.3) |

## Evo 2 (exploratory, 98.75% interval)

Evo 2 7B was scored on every case variant after the other results had been seen (PROTOCOL.md section 12). Lookups + Evo 2 minus Lookups combined, top 1, with the 95% and the Bonferroni-adjusted 98.75% interval:

| Analysis | SNV gain, 95% | SNV 98.75% | SNV p | Indel gain, 95% | Indel 98.75% | Indel p |
|---|---|---|---|---|---|---|
| main | +1.0 (-3.5 to +5.5) | (-4.5 to +6.0) | 0.760 | +1.7 (-5.9 to +9.9) | (-9.0 to +12.8) | 0.746 |
| bg1pct | +0.5 (-3.5 to +5.5) | (-4.5 to +7.0) | 0.970 | +7.5 (0.0 to +15.8) | (-2.4 to +18.1) | 0.064 |
| nomoi | +5.5 (+1.0 to +10.5) | (+0.1 to +12.0) | 0.014 | +5.8 (0.0 to +12.1) | (-0.9 to +13.8) | 0.058 |
| all1kg | -1.0 (-5.5 to +3.5) | (-7.0 to +4.5) | 0.758 | +5.0 (-2.5 to +12.1) | (-4.1 to +14.2) | 0.224 |

- SNVs: no gain under the main filters. The one 98.75% interval that excludes zero, nomoi, is discussed above; on the 500 the same check gives 0.0.
- Indels: no top-1 gain under any filter at 98.75%. At top 5, adding Evo 2 gains +4.2 (+0.8 to +8.2), p = 0.012, 98.75% (+0.8 to +9.5) in main, taking Lookups combined from 95.8% (91.8 to 99.2) to 100.0%; the gain is under 5 points, and top 5 is a descriptive outcome. At top 5 the 98.75% interval also excludes zero in nomoi (+5.0, +0.8 to +10.2) and all1kg (+4.2, +0.8 to +9.3), but not in bg1pct (+5.0, 0.0 to +11.6).
- Evo 2 alone ranks indels well: 81.7% (73.7 to 89.2) top 1 and 100.0% top 5 in main, against our SpliceAI run at 58.3% (49.6 to 66.4) and 90.8% (84.5 to 96.0). On SNVs it reaches 51.5% (43.5 to 59.5), about the level of CADD (51.0%, 43.0 to 59.5).

Against the first run (`../../rare-disease/README.md`, `REPLICATION.md`): on its 120 indel cases Evo 2 added +15.8 (+5.0 to +26.8) at top 1 with no frequency filter (median 37 background indels, lookups 61.7%) and +8.3 (+0.9 to +16.3) with a gnomAD AF filter below 0.1% (median 15, lookups 71.7%). Here, with grpmax and inheritance filters, a median list of 4 indels and lookups at 84.2%, the gain is +1.7 (-5.9 to +9.9). On SNVs the first run found +3.5 (-1.5 to +9.0), and +1.5 (-3.5 to +7.0) with the AF filter. The indel cases are a new draw under a different design, so the two runs are not paired.

## Per group, top 1 (main analysis)

### SNV

Groups with fewer than 20 cases are marked small.

| Group | Cases | Lookups combined | Lookups + AVI | Lookups + Evo 2 | Rule-based | Rule + AVI | AVI | Evo 2 |
|---|---|---|---|---|---|---|---|---|
| Canonical splice | 40 | 82.5 (68.3 to 95.6) | 85.0 (70.5 to 97.5) | 87.5 (73.2 to 97.6) | 97.5 (91.7 to 100.0) | 97.5 (91.7 to 100.0) | 87.5 (71.1 to 100.0) | 70.0 (52.4 to 86.8) |
| Missense | 60 | 85.0 (73.3 to 93.7) | 83.3 (70.8 to 93.5) | 83.3 (71.7 to 93.2) | 85.0 (75.4 to 93.2) | 86.7 (78.0 to 94.5) | 78.3 (66.1 to 89.5) | 43.3 (29.8 to 59.0) |
| Non-coding | 70 | 77.1 (64.0 to 88.1) | 88.6 (78.1 to 97.0) | 78.6 (67.7 to 88.2) | 60.0 (44.2 to 75.8) | 71.4 (58.2 to 84.1) | 67.1 (55.1 to 78.1) | 40.0 (28.6 to 50.8) |
| Non-coding: intronic | 38 | 89.5 (76.3 to 100.0) | 97.4 (91.7 to 100.0) | 86.8 (73.0 to 97.4) | 63.2 (44.7 to 82.5) | 79.0 (61.9 to 94.3) | 68.4 (53.7 to 82.4) | 36.8 (20.0 to 52.8) |
| Non-coding: splice region | 23 | 69.6 (47.6 to 90.5) | 87.0 (69.6 to 100.0) | 82.6 (66.6 to 95.8) | 60.9 (40.9 to 81.0) | 69.6 (47.6 to 90.5) | 65.2 (43.5 to 86.4) | 52.2 (34.8 to 69.6) |
| Non-coding: synonymous | 6 (small) | 66.7 (33.3 to 100.0) | 83.3 (50.0 to 100.0) | 50.0 (16.7 to 83.3) | 66.7 (16.7 to 100.0) | 66.7 (33.3 to 100.0) | 83.3 (50.0 to 100.0) | 33.3 (0.0 to 66.7) |
| Non-coding: UTR | 3 (small) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 33.3 (0.0 to 50.0) | 0.0 (0.0 to 0.0) |
| Nonsense | 30 | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 96.7 (90.0 to 100.0) | 70.0 (53.3 to 85.7) |

Single scores:

| Group | Cases | CADD | SpliceAI | REVEL | AlphaMissense | AVI | GPN-MSA | GPN-Star-V | GPN-Star-M | Evo 2 |
|---|---|---|---|---|---|---|---|---|---|---|
| Canonical splice | 40 | 77.5 (64.1 to 88.6) | 90.0 (80.0 to 97.6) | 2.5 (0.0 to 8.3) | 0.0 (0.0 to 0.0) | 87.5 (71.1 to 100.0) | 55.0 (37.5 to 73.2) | 77.5 (60.5 to 92.3) | 77.5 (64.1 to 89.8) | 70.0 (52.4 to 86.8) |
| Missense | 60 | 61.7 (48.3 to 74.6) | 15.0 (6.9 to 24.1) | 90.0 (81.8 to 96.6) | 78.3 (67.2 to 88.9) | 78.3 (66.1 to 89.5) | 60.0 (46.8 to 74.1) | 55.0 (41.1 to 69.7) | 58.3 (45.6 to 71.2) | 43.3 (29.8 to 59.0) |
| Non-coding | 70 | 7.1 (1.4 to 13.9) | 81.4 (68.9 to 92.5) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 67.1 (55.1 to 78.1) | 30.0 (17.3 to 42.4) | 24.3 (12.7 to 35.6) | 17.1 (9.4 to 25.4) | 40.0 (28.6 to 50.8) |
| Non-coding: intronic | 38 | 10.5 (2.4 to 22.2) | 86.8 (71.7 to 100.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 68.4 (53.7 to 82.4) | 36.8 (21.9 to 53.0) | 28.9 (15.6 to 43.2) | 23.7 (10.8 to 37.8) | 36.8 (20.0 to 52.8) |
| Non-coding: splice region | 23 | 0.0 (0.0 to 0.0) | 78.3 (60.0 to 95.2) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 65.2 (43.5 to 86.4) | 30.4 (10.0 to 50.0) | 26.1 (8.7 to 45.8) | 8.7 (0.0 to 20.8) | 52.2 (34.8 to 69.6) |
| Non-coding: synonymous | 6 (small) | 16.7 (0.0 to 50.0) | 100.0 (100.0 to 100.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 83.3 (50.0 to 100.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 33.3 (0.0 to 66.7) |
| Non-coding: UTR | 3 (small) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 33.3 (0.0 to 50.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 33.3 (0.0 to 50.0) | 0.0 (0.0 to 0.0) |
| Nonsense | 30 | 96.7 (89.3 to 100.0) | 16.7 (5.7 to 31.2) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 96.7 (90.0 to 100.0) | 40.0 (24.1 to 56.3) | 80.0 (63.3 to 96.5) | 76.7 (60.6 to 90.6) | 70.0 (53.3 to 85.7) |

Lookups + GPN scores:

| Group | Cases | Lookups + GPN-MSA | Lookups + GPN-Star-V | Lookups + GPN-Star-M |
|---|---|---|---|---|
| Canonical splice | 40 | 85.0 (70.5 to 97.5) | 82.5 (65.0 to 97.5) | 85.0 (71.4 to 95.5) |
| Missense | 60 | 83.3 (72.9 to 92.9) | 85.0 (74.6 to 94.5) | 85.0 (73.8 to 95.0) |
| Non-coding | 70 | 72.9 (62.0 to 83.3) | 77.1 (65.8 to 88.2) | 78.6 (67.6 to 88.9) |
| Non-coding: intronic | 38 | 86.8 (76.9 to 95.1) | 92.1 (83.3 to 100.0) | 89.5 (79.5 to 97.4) |
| Non-coding: splice region | 23 | 65.2 (45.5 to 82.6) | 65.2 (45.8 to 83.3) | 78.3 (58.3 to 95.5) |
| Non-coding: synonymous | 6 (small) | 50.0 (16.7 to 83.3) | 66.7 (33.3 to 100.0) | 50.0 (16.7 to 83.3) |
| Non-coding: UTR | 3 (small) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) |
| Nonsense | 30 | 96.7 (89.3 to 100.0) | 96.7 (89.3 to 100.0) | 96.7 (89.3 to 100.0) |

Primary and key secondary by group (top 1 gain, 95% paired interval):

| Group | Cases | Lookups + AVI minus Lookups combined | Rule + AVI minus Rule-based | Lookups combined minus Rule-based | Lookups + Evo 2 minus Lookups combined (exploratory) |
|---|---|---|---|---|---|
| Canonical splice | 40 | +2.5 (-12.5 to +17.1), p = 0.908 | 0.0 (-7.1 to +7.5), p = 1.000 | -15.0 (-29.3 to -2.5), p = 0.034 | +5.0 (-5.3 to +17.5), p = 0.560 |
| Missense | 60 | -1.7 (-10.8 to +6.6), p = 0.866 | +1.7 (0.0 to +5.3), p = 0.752 | 0.0 (-8.3 to +7.8), p = 1.000 | -1.7 (-6.8 to +3.6), p = 0.792 |
| Non-coding | 70 | +11.4 (+2.7 to +22.4), p = 0.008 | +11.4 (+2.5 to +21.6), p = 0.030 | +17.1 (+6.2 to +28.3), p = 0.002 | +1.4 (-6.9 to +10.1), p = 0.836 |
| Non-coding: intronic | 38 | +7.9 (-2.9 to +22.0), p = 0.332 | +15.8 (+4.9 to +29.7), p = 0.008 | +26.3 (+11.6 to +42.9), p = 0.001 | -2.6 (-13.2 to +5.4), p = 0.754 |
| Non-coding: splice region | 23 | +17.4 (+4.0 to +36.4), p = 0.034 | +8.7 (0.0 to +20.0), p = 0.218 | +8.7 (0.0 to +20.0), p = 0.218 | +13.0 (-4.3 to +31.8), p = 0.246 |
| Non-coding: synonymous | 6 | +16.7 (0.0 to +50.0), p = 0.706 | 0.0 (-50.0 to +50.0), p = 1.000 | 0.0 (-50.0 to +50.0), p = 1.000 | -16.7 (-50.0 to 0.0), p = 0.676 |
| Non-coding: UTR | 3 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 |
| Nonsense | 30 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 | 0.0 (0.0 to 0.0), p = 1.000 |

Both AVI gains come from the non-coding group (primary +11.4 (+2.7 to +22.4), p = 0.008, key secondary +11.4 (+2.5 to +21.6), p = 0.030). Evo 2 gains nothing in any group; its largest estimate, splice region +13.0 (-4.3 to +31.8), p = 0.246, has an interval wide enough to include zero. Neither gains on nonsense; on canonical splice and missense the gains are near 0 with wide intervals.

### Indel

| Group | Cases | SpliceAI | Evo 2 | Rule-based | Lookups combined | Lookups + Evo 2 |
|---|---|---|---|---|---|---|
| Canonical splice | 30 | 100.0 (100.0 to 100.0) | 93.3 (82.7 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 96.7 (89.3 to 100.0) |
| Frameshift | 40 | 32.5 (15.8 to 50.0) | 87.5 (75.0 to 97.6) | 87.5 (74.4 to 97.8) | 92.5 (82.9 to 100.0) | 90.0 (78.0 to 100.0) |
| Inframe | 30 | 46.7 (27.6 to 65.5) | 83.3 (66.7 to 96.9) | 70.0 (50.0 to 88.9) | 76.7 (59.4 to 92.9) | 86.7 (74.1 to 96.9) |
| Non-coding | 20 | 65.0 (42.9 to 85.0) | 50.0 (30.0 to 68.4) | 30.0 (10.5 to 52.6) | 55.0 (35.0 to 73.7) | 60.0 (40.0 to 79.0) |
| Non-coding: intronic | 13 (small) | 92.3 (75.0 to 100.0) | 69.2 (43.7 to 91.7) | 38.5 (13.3 to 66.7) | 76.9 (56.2 to 100.0) | 76.9 (53.8 to 100.0) |
| Non-coding: UTR | 7 (small) | 14.3 (0.0 to 42.9) | 14.3 (0.0 to 42.9) | 14.3 (0.0 to 42.9) | 14.3 (0.0 to 42.9) | 28.6 (0.0 to 71.4) |

| Group | Cases | Lookups combined minus Rule-based | Lookups + Evo 2 minus Lookups combined, top 1 | Lookups + Evo 2 minus Lookups combined, top 5 |
|---|---|---|---|---|
| Canonical splice | 30 | 0.0 (0.0 to 0.0), p = 1.000 | -3.3 (-10.7 to 0.0), p = 0.758 | 0.0 (0.0 to 0.0), p = 1.000 |
| Frameshift | 40 | +5.0 (0.0 to +16.7), p = 0.756 | -2.5 (-16.7 to +9.5), p = 0.900 | 0.0 (0.0 to 0.0), p = 1.000 |
| Inframe | 30 | +6.7 (0.0 to +16.1), p = 0.252 | +10.0 (-3.9 to +26.7), p = 0.344 | 0.0 (0.0 to 0.0), p = 1.000 |
| Non-coding | 20 | +25.0 (+5.6 to +42.1), p = 0.010 | +5.0 (-10.5 to +25.0), p = 0.784 | +25.0 (+5.6 to +44.5), p = 0.002 |
| Non-coding: intronic | 13 | +38.5 (+15.4 to +61.5), p = 0.004 | 0.0 (-23.1 to +21.5), p = 1.000 | +7.7 (0.0 to +25.0), p = 0.682 |
| Non-coding: UTR | 7 | 0.0 (0.0 to 0.0), p = 1.000 | +14.3 (0.0 to +42.9), p = 0.610 | +57.1 (+28.6 to +100.0), p = 0.006 |

Intronic (13) and UTR (7) are small. The rule gives intronic indels only a SpliceAI level, so it ranks them below Lookups combined (38.5%, 13.3 to 66.7, against 76.9%, 56.2 to 100.0), and below SpliceAI alone (92.3%, 75.0 to 100.0).

Evo 2 alone ranks frameshift (87.5%, 75.0 to 97.6) and inframe (83.3%, 66.7 to 96.9) diagnostics first about as often as Lookups combined, and is weaker on non-coding (50.0%, 30.0 to 68.4). Added to the lookups it gains at top 1 only on inframe, +10.0 (-3.9 to +26.7), p = 0.344. Its whole top-5 gain (+4.2, five cases) is in the 20 non-coding cases: +25.0 (+5.6 to +44.5), p = 0.002.

## Diagnostic found in gnomAD or absent (section 9)

Found means the diagnostic is in full gnomAD v3.1.2 genomes with AC above 0 (`data/gnomad_full.tsv`, then `data/pool/gnomad_full.tsv`). SNV: 52 found, 148 absent. Indel: 29 found, 91 absent. Every diagnostic has a status. The contrasts for the two strata are given only in the `gnomad_split` block, where both use one shared draw of people (seed label "<study>/<analysis>/gnomad split"), so absent minus found is paired.

| SNV analysis | Contrast | Found | Absent | Absent minus found | Flagged |
|---|---|---|---|---|---|
| main | primary | +5.8 (-2.1 to +15.1), p = 0.236 | +3.4 (-2.8 to +10.0), p = 0.348 | -2.4 (-12.2 to +6.0), p = 0.598 | No |
| main | key secondary | 0.0 (-6.7 to +5.0), p = 1.000 | +6.1 (+1.9 to +11.5), p = 0.010 | +6.1 (-0.6 to +13.8), p = 0.082 | No |
| bg1pct | primary | +3.9 (-13.1 to +18.8), p = 0.680 | 0.0 (-9.4 to +9.3), p = 1.000 | -3.9 (-19.0 to +13.6), p = 0.596 | No |
| bg1pct | key secondary | 0.0 (-4.7 to +5.3), p = 1.000 | +5.4 (+0.7 to +10.4), p = 0.026 | +5.4 (-1.2 to +12.7), p = 0.120 | No |
| nomoi | primary | +9.6 (0.0 to +20.4), p = 0.094 | +4.0 (-4.7 to +13.2), p = 0.398 | -5.6 (-17.1 to +6.4), p = 0.392 | No |
| nomoi | key secondary | +3.9 (-3.9 to +11.5), p = 0.430 | +8.1 (+2.9 to +13.8), p = 0.002 | +4.3 (-4.1 to +13.3), p = 0.316 | No |
| all1kg | primary | +7.7 (-2.0 to +18.6), p = 0.180 | +7.4 (+0.7 to +15.7), p = 0.040 | -0.3 (-12.5 to +11.9), p = 0.934 | No |
| all1kg | key secondary | -1.9 (-6.7 to 0.0), p = 0.714 | +6.8 (+2.6 to +12.3), p = 0.001 | +8.7 (+3.2 to +15.5), p = 0.001 | Yes |

- **Primary:** not flagged in any analysis. The gain from AVI in Lookups + AVI is no larger in absent cases; in main it is larger in found cases.
- **Key secondary:** Rule + AVI gains nothing in found cases in main, bg1pct and all1kg, and gains with an interval above 0 in absent cases in every analysis. `case_benchmark.py` flags it in all1kg (absent minus found +8.7 (+3.2 to +15.5)), and in main the difference nearly excludes 0 (+6.1 (-0.6 to +13.8)). On the 200, the all1kg check flags the key secondary gain as possibly due to gnomAD presence, not to what a lab would see; main, bg1pct and nomoi do not flag it. On the 500 no analysis flags it (see `analysis_replication.md`).
- Top 1 in main, found against absent: AVI 73.1 (60.3 to 86.0) against 81.1 (71.6 to 89.8); CADD 51.9 (40.0 to 63.8) against 50.7 (41.6 to 60.1); Rule-based 78.8 (67.3 to 90.6) against 81.8 (73.9 to 89.3); Lookups combined 80.8 (68.0 to 90.9) against 85.1 (76.4 to 93.3).
- Indel: no AVI contrast, so the split has cases only. Evo 2 alone, found against absent: 82.8 (66.7 to 96.3) and 81.3 (72.6 to 90.0). Top 1 in main, found against absent: SpliceAI 65.5 (50.0 to 80.0) and 56.0 (44.9 to 66.7), Rule-based 82.8 (66.6 to 96.7) and 74.7 (63.2 to 85.7), Lookups combined 86.2 (72.4 to 96.7) and 83.5 (74.4 to 91.8).

## AVI training chromosomes (snv heldout)

Cases split by whether the diagnostic is on an AVI training chromosome (1, 4, 7, 8, 10, 13, 15); each ranked only among background variants in the same split. Both splits use the same draws of people.

| Split | Cases | Median list | AVI | Lookups combined | Lookups + AVI | Rule-based | Rule + AVI | Evo 2 | Lookups + Evo 2 | Primary gain | Key secondary gain |
|---|---|---|---|---|---|---|---|---|---|---|---|
| training | 38 | 7 | 81.6 (67.6 to 93.1) | 94.7 (87.8 to 100.0) | 97.4 (91.2 to 100.0) | 89.5 (80.0 to 97.4) | 92.1 (84.2 to 100.0) | 76.3 (60.5 to 89.3) | 94.7 (88.1 to 100.0) | +2.6 (-6.2 to +10.8), p = 0.768 | +2.6 (0.0 to +8.9), p = 0.738 |
| other | 162 | 19 | 84.0 (76.1 to 90.8) | 88.3 (80.8 to 94.5) | 91.4 (84.7 to 96.7) | 82.7 (74.7 to 90.5) | 87.0 (79.5 to 93.2) | 54.3 (46.1 to 62.7) | 88.3 (81.8 to 94.2) | +3.1 (-0.6 to +7.5), p = 0.182 | +4.3 (+0.6 to +8.2), p = 0.022 |

Difference, training minus other: primary -0.4 (-9.8 to +8.6), p = 0.892; key secondary -1.7 (-7.5 to +5.6), p = 0.604. `possibly_inflated_by_training_overlap` is false. Evo 2 alone is higher on the training split (76.3%, 60.5 to 89.3, against 54.3%, 46.1 to 62.7), where lists are shorter; the split is AVI's training chromosomes and says nothing about Evo 2. The training split has 38 cases with a median list of 7, so one or two cases move these numbers by several points.

## EVEE (snv evee, descriptive only)

Agreement with ClinVar labels on in-sample variants, not diagnostic performance (section 10). bg1pct filters, lists restricted to variants with an EVEE row (7,084 of 19,561). Every one of the 200 diagnostics has a row, so no case is dropped; median list 29.5. Not to be compared with the main analysis.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 44.5 (37.5 to 51.5) | 74.0 (68.0 to 80.5) | 91.0 (86.0 to 95.5) |
| SpliceAI | 50.0 (43.5 to 57.0) | 63.5 (58.0 to 68.5) | 66.0 (61.0 to 71.0) |
| REVEL | 26.5 (21.5 to 31.5) | 31.5 (26.0 to 37.0) | 33.0 (27.5 to 38.5) |
| AlphaMissense | 24.0 (19.5 to 28.5) | 28.0 (23.5 to 32.5) | 28.5 (24.0 to 33.0) |
| AVI | 74.5 (66.5 to 81.5) | 93.5 (89.5 to 96.5) | 97.5 (95.0 to 99.5) |
| GPN-MSA | 41.5 (33.5 to 49.0) | 74.5 (69.5 to 79.0) | 87.0 (83.0 to 91.0) |
| GPN-Star-V | 51.5 (44.0 to 58.5) | 82.0 (77.5 to 87.0) | 92.0 (88.0 to 95.5) |
| GPN-Star-M | 42.5 (34.0 to 51.0) | 81.0 (76.0 to 85.5) | 92.0 (88.5 to 95.5) |
| Evo 2 | 47.0 (39.0 to 55.0) | 80.0 (73.5 to 86.0) | 91.0 (87.0 to 94.5) |
| EVEE | 81.5 (75.5 to 88.5) | 97.5 (95.0 to 99.5) | 98.5 (96.5 to 100.0) |
| Lookups combined | 82.0 (74.5 to 88.5) | 96.5 (93.0 to 99.0) | 97.5 (94.0 to 99.5) |

EVEE and Lookups combined are close (81.5%, 75.5 to 88.5; 82.0%, 74.5 to 88.5), as on the 500. Evo 2 alone reaches 47.0% (39.0 to 55.0).

## Checks

- **Invariance.** `tests/test_invariance.py` reruns an analysis with one score dropped and compares every row that does not involve it (`logs/cb_invariance.log`): snv main without GPN-Star-M, 1,708 rows compared, 0 problems; without AVI, 1,452 rows, 0 problems; snv heldout without GPN-Star-M, 278 rows, 0 problems. Indel main without GPN-Star-M: 348 rows compared, 0 problems (GPN-Star-M has no indel rows). That log predates Evo 2. Adding Evo 2 was checked against tag `amp-2026-v5` instead: every value in the 16 result files that does not involve Evo 2 is identical.
- **Lookups combined against the best single tool (SNV, top 1).** Point estimates, with overlapping intervals in every analysis: main 84.0% (75.5 to 91.5) against AVI 79.0% (71.0 to 86.0), higher; bg1pct 77.0% (68.0 to 85.0) against AVI 70.0% (61.5 to 78.0), higher; nomoi 68.5% (57.5 to 79.0) against AVI 66.5% (58.5 to 74.0), higher; all1kg 77.0% (68.0 to 85.5) against AVI 77.5% (69.0 to 85.5), lower; evee 82.0% (74.5 to 88.5) against EVEE 81.5% (75.5 to 88.5). In the indel study it beats SpliceAI in every analysis, and is level with or above Evo 2 alone: main 84.2% (75.4 to 92.4) against 81.7% (73.7 to 89.2); bg1pct 75.0% (65.5 to 83.6) against 75.0% (67.5 to 81.7); nomoi 81.7% (73.0 to 89.9) against 78.3% (70.8 to 86.2); all1kg 83.3% (75.4 to 90.8) against 81.7% (73.6 to 89.1). Within groups it falls below a single tool in places: SNV canonical splice 82.5% (68.3 to 95.6) against SpliceAI 90.0% (80.0 to 97.6), and indel intronic 76.9% (56.2 to 100.0) against SpliceAI 92.3% (75.0 to 100.0) (13 cases).
- **Charts.** `figures/rd_charts.py` draws gains_primary, top1_by_group, gains_added, gains_gnomad, avi_first_vs_second and evo2_gains from these files and the replication files, into `figures/out/` (not committed); with `--readme` it also copies the README's charts into `figures/` (committed).
