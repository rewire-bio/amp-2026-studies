# Results on the 200 SNV and 120 indel cases (cross-validated)

These are the 200 SNV and 120 indel cases, with every model fitted by 5-fold cross-validation grouped by person. They are descriptive. They are not the primary outcome, which is the locked models on the 500 replication cases (PROTOCOL.md section 8). No verdict is read off these numbers: `case_benchmark.py` applies the reading rule only to the locked primary on the 500.

Made by `case_benchmark.py` on the real case and score files, 2 November 2026, after the code fixes (pairwise rule ties, type drop after the inheritance filter, gnomAD split) and the data rebuild (indels in every background, multi-base indels left out of the pool). Result files: `results/snv_<analysis>.json` for main, bg1pct, nomoi, all1kg, heldout and evee, and `results/indel_<analysis>.json` for main, bg1pct, nomoi and all1kg. Every run exited 0. Log: `logs/analysis_snv_indel.log`. Intervals are 95% bootstrap intervals over people (1,000 resamples); every gain is a paired interval. Top k is the share of cases with fewer than k competitors. For single scores and models a competitor is a background variant scoring equal to or higher than the diagnostic. For the rule scores it is a background variant whose pairwise comparison with the diagnostic is TIED or BEATS (PROTOCOL.md section 7).

In the indel study, "SpliceAI" is the score we ran, and each indel diagnostic is ranked among the person's background indels only. CADD and precomputed SpliceAI are set aside for indels by design (section 6), and REVEL, AlphaMissense, AVI and the GPN scores have no indel rows. So the indel study has SpliceAI, Rule-based and Lookups combined only, and no AVI contrast.

## Headline (main analysis)

- SNV, 200 cases, median list 26 variants (diagnostic included). Top 1: Lookups combined 0.840, Rule-based 0.810, AVI alone 0.790, Lookups + AVI 0.880, Rule + AVI 0.855.
- Primary contrast, cross-validated: Lookups + AVI minus Lookups combined, top 1 +0.040 [-0.020, 0.095].
- Key secondary: Rule + AVI minus Rule-based, top 1 +0.045 [0.005, 0.080].
- GPN family added to Lookups combined: GPN-MSA -0.020 [-0.070, 0.025], GPN-Star-V -0.005 [-0.055, 0.040], GPN-Star-M +0.005 [-0.040, 0.045].
- The key secondary gain sits in the cases whose diagnostic is absent from gnomAD (+0.061 [0.020, 0.116]; 0.000 in found cases). See the gnomAD split below.
- Indel, 120 cases, median list 4 indels. Top 1: Lookups combined 0.842, Rule-based 0.767, SpliceAI 0.583. Lookups combined minus Rule-based +0.075 [0.033, 0.126].

## What changed against the previous numbers, and why

Three runs of snv main: the previous file (old code, old case files), the new code on the old case files, and the new code on the new case files. The middle run is not one of the result files: it is `case_benchmark.py` as it is now, with RDV2_DATA pointed at a copy of data/ whose snv and indel folders are the previous case and score files kept by the data rebuild, and RDV2_RESULTS at a scratch folder. The old snv case files are the new ones without the 4,524 background indel rows; the 200 diagnostics and every background SNV row are identical.

| Top 1, snv main | Previous | New code, old case files | New code, new case files |
|---|---|---|---|
| Lookups combined | 0.810 | 0.810 | 0.840 |
| Rule-based | 0.825 | 0.810 | 0.810 |
| Lookups + AVI | 0.865 | 0.865 | 0.880 |
| Rule + AVI | 0.860 | 0.855 | 0.855 |
| Primary (Lookups + AVI minus Lookups combined) | +0.055 [-0.015, 0.115] | +0.055 [-0.015, 0.115] | +0.040 [-0.020, 0.095] |
| Key secondary (Rule + AVI minus Rule-based) | +0.035 [-0.005, 0.070] | +0.045 [0.005, 0.080] | +0.045 [0.005, 0.080] |
| Variants in all lists | 5,369 | 5,369 | 5,401 |

- **Rule-based 0.825 to 0.810, Rule + AVI 0.860 to 0.855, key secondary +0.035 to +0.045: the pairwise tie rule.** These move with the code change alone and not with the data. Scores from different tools are no longer compared; CADD decides between them. Rule-based top 5 also fell, from 0.980 to 0.960. Rule-based lost more top 1 than Rule + AVI did, so the key secondary grew and its interval now excludes 0. Rule-based is also lower in nomoi (0.695 to 0.665), where the lists did not change, so that is the rule change alone. In bg1pct (0.745 to 0.715) and all1kg (0.820 to 0.795) the lists also grew, so both fixes act there. By group in main it is the non-coding diagnostics: Rule-based 0.66 to 0.60, intronic 0.74 to 0.63.
- **Lookups combined 0.810 to 0.840, Lookups + AVI 0.865 to 0.880, primary +0.055 to +0.040: the inheritance filter fix.** Counting each person's indels in the two-variant rule for biallelic genes keeps 32 more background SNVs in 28 cases (788 background indels enter the count; 6,189 rows after the filter, 5,401 after the indels are dropped). Only 2 cases have an added SNV that scores at or above the diagnostic under Lookups combined, and neither lost top 1. The gain comes from the refit instead: `GroupKFold` assigns people to folds by row count, so the 32 rows moved 25 of the 50 people to a different fold. Lookups combined gained 6 cases and lost 0; Lookups + AVI gained 4 and lost 1. A 3-point move from fold assignment alone is inside every interval here, and is one reason these cross-validated numbers are descriptive. The nomoi run confirms the cause: it has no inheritance filter, its lists are unchanged (9,816 variants), and Lookups combined (0.685) and Lookups + AVI (0.740) are unchanged there too.
- **The locked-fit fix** (training set filtered under main, or bg1pct for bg1pct and evee) applies to the replication study only. It changes nothing here.
- **The gnomAD split** is new; it adds strata and a block, and changes no existing row.
- **Model classes** are unchanged: synonymous and UTR merge into other for SNVs; UTR merges into other for indels.

Indel main is a new sample, not a before and after. The rebuild left out multi-base indels from the pool, and 119 of the 120 cases now have a different diagnostic (same person per case). The new code on the old indel case files isolates the rule change: Rule-based 0.733 to 0.717, tie share 0.117 to 0.083, Lookups combined unchanged at 0.758. The rest is the new cases: Rule-based 0.717 to 0.767, Lookups combined 0.758 to 0.842, SpliceAI 0.583 in both. Lookups combined minus Rule-based went from +0.025 [-0.017, 0.076] to +0.075 [0.033, 0.126]. Indel diagnostics without a SpliceAI score fell from 17 to 9 (6 of the 7 UTR diagnostics, and 1 each of canonical splice, frameshift and intronic).

## Cases and list lengths

Median list length per case, diagnostic included. No case lost its diagnostic under any filter (`cases_lost_diagnostic` is 0 in every run).

| Analysis | SNV cases | SNV median list | SNV variants in lists | Indel cases | Indel median list | Indel variants in lists |
|---|---|---|---|---|---|---|
| main | 200 | 26 | 5,401 | 120 | 4 | 576 |
| bg1pct | 200 | 78.5 | 19,561 | 120 | 14 | 1,955 |
| nomoi | 200 | 51.5 | 9,816 | 120 | 7 | 781 |
| all1kg | 200 | 25.5 | 5,581 | 120 | 4 | 610 |
| heldout | 200 | 26 | 5,401 | not run (SNV only) | | |
| evee | 200 | 29.5 | 7,084 | not run (SNV only) | | |

- The inheritance filter now counts both variant types, then keeps only the study's own type. SNV main: 10,956 rows after the frequency filter, 6,189 after the inheritance filter (788 of them background indels), 5,401 SNVs ranked. Indel main: 6,410, then 3,615 (3,039 background SNVs), 576 indels ranked.
- The nomoi medians match the build reports, which give list lengths after the frequency filter only: SNV 51.5 is 50.5 background variants (build_snv.md 50.5), indel 7 is 6 background indels (build_indel.md 6).

## Model class merges

- SNV (from the 200): canonical splice 40, intronic 38, missense 60, nonsense 30, splice region 23, synonymous 6, UTR 3. Merged into other: synonymous, UTR. Model classes: canonical splice, intronic, missense, nonsense, splice region, other.
- Indel (from the 120): canonical splice 30, frameshift 40, inframe 30, intronic 13, UTR 7. Merged into other: UTR. Intronic has 13 diagnostics, so it stays.
- Each list now holds only its own type's classes. The same merges hold in every analysis. `diagnostic_class_differs_from_cases_tsv` and `case_variants_without_mane_cadd_row` are 0 in every run.

## Rule tie share

Share of diagnostics with at least one background variant whose pairwise outcome is TIED.

- SNV: 0 for Rule-based and Rule + AVI in main, bg1pct, all1kg and every group of main. nomoi: Rule-based 0.010 (2 of 200), Rule + AVI 0.
- Indel: main 0.042 (5 of 120), bg1pct 0.092, nomoi 0.050, all1kg 0.042. By group in main: frameshift 0.050 (2 of 40), inframe 0.100 (3 of 30), 0 elsewhere.
- Indel ties are indels with the same points and the same SpliceAI score (often 0) as the diagnostic; indels have no CADD to break them. The previous indel tie share of 0.117 is not comparable: it was on the old diagnostics and under the old rule.

## Coverage per score (main analysis)

Variants with a score over all variants in the lists (diagnostics included), and diagnostics with a score. CADD shows nothing for indels because it is set aside for them.

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
| EVEE (evee analysis) | 7,084 (the EVEE-covered bg1pct lists, out of 19,561) | 200 | not used | not used |

The 32 added SNVs account for the small rises against the previous table (5,369 variants, SpliceAI 4,461, REVEL 1,240, AlphaMissense 1,021).

## Results by analysis, all cases

### SNV study, 200 cases

#### main (grpmax below 0.1%, inheritance filter)

200 cases, median list 26. Rule tie share: Rule-based 0.000, Rule + AVI 0.000.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 0.510 [0.430, 0.595] | 0.820 [0.760, 0.875] | 0.930 [0.900, 0.960] |
| SpliceAI | 0.535 [0.460, 0.595] | 0.670 [0.610, 0.720] | 0.705 [0.640, 0.765] |
| REVEL | 0.275 [0.230, 0.325] | 0.325 [0.275, 0.375] | 0.380 [0.315, 0.455] |
| AlphaMissense | 0.235 [0.190, 0.280] | 0.285 [0.245, 0.335] | 0.325 [0.265, 0.395] |
| AVI | 0.790 [0.710, 0.860] | 0.950 [0.920, 0.980] | 0.990 [0.975, 1.000] |
| GPN-MSA | 0.455 [0.360, 0.550] | 0.765 [0.705, 0.825] | 0.885 [0.835, 0.930] |
| GPN-Star-V | 0.525 [0.440, 0.605] | 0.855 [0.800, 0.900] | 0.920 [0.880, 0.960] |
| GPN-Star-M | 0.505 [0.425, 0.580] | 0.805 [0.740, 0.865] | 0.935 [0.895, 0.965] |
| Rule-based | 0.810 [0.735, 0.880] | 0.960 [0.925, 0.990] | 0.975 [0.945, 0.995] |
| Rule + AVI | 0.855 [0.790, 0.910] | 0.965 [0.935, 0.990] | 0.975 [0.945, 0.995] |
| Lookups combined | 0.840 [0.755, 0.915] | 0.970 [0.945, 0.995] | 0.985 [0.960, 1.000] |
| Lookups + AVI | 0.880 [0.810, 0.940] | 0.980 [0.955, 1.000] | 0.980 [0.955, 1.000] |
| Lookups + GPN-MSA | 0.820 [0.740, 0.890] | 0.960 [0.930, 0.985] | 0.970 [0.940, 0.995] |
| Lookups + GPN-Star-V | 0.835 [0.755, 0.910] | 0.970 [0.945, 0.995] | 0.975 [0.950, 0.995] |
| Lookups + GPN-Star-M | 0.845 [0.755, 0.915] | 0.970 [0.940, 0.995] | 0.980 [0.955, 1.000] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast | +0.040 [-0.020, 0.095] | +0.010 [-0.010, 0.025] | -0.005 [-0.015, 0.000] |
| Lookups combined minus Rule-based | descriptive | +0.030 [-0.045, 0.095] | +0.010 [-0.015, 0.040] | +0.010 [0.000, 0.025] |
| Lookups + AVI minus Rule-based | descriptive | +0.070 [0.000, 0.130] | +0.020 [-0.005, 0.050] | +0.005 [-0.010, 0.020] |
| Lookups + GPN-MSA minus Lookups combined | secondary | -0.020 [-0.070, 0.025] | -0.010 [-0.025, 0.000] | -0.015 [-0.035, 0.000] |
| Lookups + GPN-Star-V minus Lookups combined | secondary | -0.005 [-0.055, 0.040] | +0.000 [-0.015, 0.015] | -0.010 [-0.025, 0.000] |
| Lookups + GPN-Star-M minus Lookups combined | secondary | +0.005 [-0.040, 0.045] | +0.000 [-0.015, 0.015] | -0.005 [-0.015, 0.000] |
| Rule + AVI minus Rule-based | key secondary | +0.045 [0.005, 0.080] | +0.005 [-0.010, 0.025] | +0.000 [0.000, 0.000] |

#### bg1pct (grpmax below 1%, inheritance filter)

200 cases, median list 78.5. Rule tie share: Rule-based 0.000, Rule + AVI 0.000.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 0.410 [0.335, 0.485] | 0.670 [0.600, 0.740] | 0.800 [0.730, 0.865] |
| SpliceAI | 0.470 [0.400, 0.545] | 0.615 [0.560, 0.670] | 0.640 [0.585, 0.695] |
| REVEL | 0.255 [0.205, 0.305] | 0.305 [0.255, 0.360] | 0.320 [0.270, 0.375] |
| AlphaMissense | 0.225 [0.180, 0.275] | 0.280 [0.235, 0.330] | 0.285 [0.240, 0.335] |
| AVI | 0.700 [0.615, 0.780] | 0.930 [0.890, 0.965] | 0.955 [0.920, 0.985] |
| GPN-MSA | 0.345 [0.275, 0.425] | 0.620 [0.555, 0.680] | 0.760 [0.720, 0.800] |
| GPN-Star-V | 0.450 [0.380, 0.520] | 0.710 [0.650, 0.770] | 0.845 [0.800, 0.890] |
| GPN-Star-M | 0.380 [0.300, 0.460] | 0.680 [0.620, 0.735] | 0.835 [0.790, 0.880] |
| Rule-based | 0.715 [0.625, 0.795] | 0.935 [0.895, 0.970] | 0.965 [0.935, 0.990] |
| Rule + AVI | 0.755 [0.680, 0.825] | 0.955 [0.920, 0.980] | 0.960 [0.930, 0.985] |
| Lookups combined | 0.770 [0.680, 0.850] | 0.955 [0.920, 0.985] | 0.965 [0.935, 0.990] |
| Lookups + AVI | 0.780 [0.695, 0.855] | 0.955 [0.920, 0.980] | 0.960 [0.925, 0.985] |
| Lookups + GPN-MSA | 0.765 [0.685, 0.840] | 0.955 [0.920, 0.980] | 0.960 [0.925, 0.985] |
| Lookups + GPN-Star-V | 0.770 [0.695, 0.840] | 0.950 [0.910, 0.980] | 0.950 [0.910, 0.980] |
| Lookups + GPN-Star-M | 0.775 [0.705, 0.840] | 0.945 [0.905, 0.975] | 0.955 [0.920, 0.980] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast | +0.010 [-0.075, 0.095] | +0.000 [-0.020, 0.020] | -0.005 [-0.030, 0.015] |
| Lookups combined minus Rule-based | descriptive | +0.055 [-0.025, 0.135] | +0.020 [-0.010, 0.055] | +0.000 [-0.020, 0.020] |
| Lookups + AVI minus Rule-based | descriptive | +0.065 [-0.025, 0.155] | +0.020 [-0.010, 0.055] | -0.005 [-0.020, 0.010] |
| Lookups + GPN-MSA minus Lookups combined | secondary | -0.005 [-0.050, 0.040] | +0.000 [-0.015, 0.015] | -0.005 [-0.015, 0.000] |
| Lookups + GPN-Star-V minus Lookups combined | secondary | +0.000 [-0.060, 0.060] | -0.005 [-0.025, 0.010] | -0.015 [-0.035, 0.000] |
| Lookups + GPN-Star-M minus Lookups combined | secondary | +0.005 [-0.045, 0.060] | -0.010 [-0.025, 0.000] | -0.010 [-0.030, 0.010] |
| Rule + AVI minus Rule-based | key secondary | +0.040 [0.005, 0.080] | +0.020 [0.000, 0.045] | -0.005 [-0.015, 0.000] |

#### nomoi (grpmax below 0.1%, no inheritance filter)

200 cases, median list 51.5. Rule tie share: Rule-based 0.010, Rule + AVI 0.000.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 0.400 [0.325, 0.470] | 0.735 [0.665, 0.800] | 0.835 [0.780, 0.890] |
| SpliceAI | 0.460 [0.390, 0.530] | 0.630 [0.575, 0.685] | 0.665 [0.610, 0.720] |
| REVEL | 0.235 [0.180, 0.290] | 0.305 [0.255, 0.360] | 0.330 [0.280, 0.385] |
| AlphaMissense | 0.200 [0.150, 0.250] | 0.280 [0.235, 0.325] | 0.280 [0.235, 0.325] |
| AVI | 0.665 [0.585, 0.740] | 0.910 [0.860, 0.950] | 0.950 [0.920, 0.980] |
| GPN-MSA | 0.355 [0.270, 0.440] | 0.675 [0.610, 0.740] | 0.820 [0.765, 0.875] |
| GPN-Star-V | 0.440 [0.360, 0.515] | 0.775 [0.715, 0.835] | 0.870 [0.825, 0.910] |
| GPN-Star-M | 0.395 [0.320, 0.470] | 0.720 [0.650, 0.785] | 0.830 [0.775, 0.885] |
| Rule-based | 0.665 [0.580, 0.750] | 0.930 [0.890, 0.965] | 0.965 [0.930, 0.990] |
| Rule + AVI | 0.735 [0.655, 0.815] | 0.940 [0.905, 0.975] | 0.960 [0.930, 0.985] |
| Lookups combined | 0.685 [0.575, 0.790] | 0.950 [0.915, 0.980] | 0.965 [0.935, 0.990] |
| Lookups + AVI | 0.740 [0.650, 0.825] | 0.945 [0.910, 0.975] | 0.965 [0.935, 0.990] |
| Lookups + GPN-MSA | 0.710 [0.615, 0.805] | 0.935 [0.895, 0.975] | 0.955 [0.925, 0.985] |
| Lookups + GPN-Star-V | 0.710 [0.620, 0.795] | 0.945 [0.910, 0.975] | 0.965 [0.940, 0.985] |
| Lookups + GPN-Star-M | 0.680 [0.590, 0.765] | 0.930 [0.890, 0.965] | 0.955 [0.925, 0.985] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast | +0.055 [-0.025, 0.135] | -0.005 [-0.020, 0.010] | +0.000 [-0.020, 0.020] |
| Lookups combined minus Rule-based | descriptive | +0.020 [-0.070, 0.125] | +0.020 [-0.005, 0.050] | +0.000 [-0.020, 0.025] |
| Lookups + AVI minus Rule-based | descriptive | +0.075 [-0.015, 0.165] | +0.015 [-0.015, 0.050] | +0.000 [-0.025, 0.030] |
| Lookups + GPN-MSA minus Lookups combined | secondary | +0.025 [-0.015, 0.070] | -0.015 [-0.035, 0.000] | -0.010 [-0.025, 0.000] |
| Lookups + GPN-Star-V minus Lookups combined | secondary | +0.025 [-0.030, 0.080] | -0.005 [-0.015, 0.000] | +0.000 [-0.015, 0.015] |
| Lookups + GPN-Star-M minus Lookups combined | secondary | -0.005 [-0.055, 0.045] | -0.020 [-0.045, 0.000] | -0.010 [-0.025, 0.000] |
| Rule + AVI minus Rule-based | key secondary | +0.070 [0.025, 0.120] | +0.010 [0.000, 0.025] | -0.005 [-0.020, 0.010] |

#### all1kg (all of 1000 Genomes left out of the background frequency, inheritance filter)

200 cases, median list 25.5. Rule tie share: Rule-based 0.000, Rule + AVI 0.000.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 0.485 [0.405, 0.565] | 0.815 [0.750, 0.875] | 0.930 [0.895, 0.960] |
| SpliceAI | 0.525 [0.455, 0.580] | 0.670 [0.610, 0.720] | 0.680 [0.620, 0.735] |
| REVEL | 0.275 [0.230, 0.330] | 0.325 [0.275, 0.380] | 0.360 [0.300, 0.435] |
| AlphaMissense | 0.235 [0.195, 0.285] | 0.285 [0.245, 0.330] | 0.310 [0.250, 0.380] |
| AVI | 0.775 [0.690, 0.855] | 0.950 [0.915, 0.980] | 0.990 [0.975, 1.000] |
| GPN-MSA | 0.435 [0.350, 0.525] | 0.765 [0.710, 0.825] | 0.885 [0.835, 0.930] |
| GPN-Star-V | 0.525 [0.445, 0.605] | 0.850 [0.800, 0.895] | 0.915 [0.870, 0.955] |
| GPN-Star-M | 0.500 [0.425, 0.570] | 0.800 [0.735, 0.855] | 0.930 [0.890, 0.965] |
| Rule-based | 0.795 [0.705, 0.865] | 0.960 [0.925, 0.990] | 0.975 [0.945, 0.995] |
| Rule + AVI | 0.840 [0.770, 0.900] | 0.965 [0.935, 0.990] | 0.975 [0.945, 0.995] |
| Lookups combined | 0.770 [0.680, 0.855] | 0.975 [0.950, 0.995] | 0.985 [0.960, 1.000] |
| Lookups + AVI | 0.845 [0.765, 0.910] | 0.975 [0.945, 0.995] | 0.975 [0.945, 0.995] |
| Lookups + GPN-MSA | 0.790 [0.705, 0.865] | 0.970 [0.940, 0.995] | 0.980 [0.955, 1.000] |
| Lookups + GPN-Star-V | 0.825 [0.745, 0.895] | 0.980 [0.955, 1.000] | 0.990 [0.975, 1.000] |
| Lookups + GPN-Star-M | 0.825 [0.750, 0.895] | 0.970 [0.940, 0.995] | 0.980 [0.960, 0.995] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups + AVI minus Lookups combined | primary contrast | +0.075 [0.020, 0.145] | +0.000 [-0.020, 0.020] | -0.010 [-0.025, 0.000] |
| Lookups combined minus Rule-based | descriptive | -0.025 [-0.110, 0.055] | +0.015 [-0.010, 0.045] | +0.010 [0.000, 0.025] |
| Lookups + AVI minus Rule-based | descriptive | +0.050 [-0.025, 0.120] | +0.015 [-0.010, 0.045] | +0.000 [-0.020, 0.020] |
| Lookups + GPN-MSA minus Lookups combined | secondary | +0.020 [-0.020, 0.060] | -0.005 [-0.015, 0.000] | -0.005 [-0.015, 0.000] |
| Lookups + GPN-Star-V minus Lookups combined | secondary | +0.055 [0.005, 0.110] | +0.005 [0.000, 0.015] | +0.005 [0.000, 0.015] |
| Lookups + GPN-Star-M minus Lookups combined | secondary | +0.055 [0.015, 0.095] | -0.005 [-0.015, 0.000] | -0.005 [-0.020, 0.010] |
| Rule + AVI minus Rule-based | key secondary | +0.045 [0.010, 0.085] | +0.005 [-0.010, 0.020] | +0.000 [0.000, 0.000] |


### Indel study, 120 cases, ranked among background indels

Rule-based top 5 is below 1 although the median list holds 4 indels: a diagnostic whose class needs SpliceAI and has none ranks below every scored indel.

#### main (grpmax below 0.1%, inheritance filter)

120 cases, median list 4. Rule tie share: Rule-based 0.042.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 0.583 [0.496, 0.664] | 0.908 [0.845, 0.960] | 0.992 [0.974, 1.000] |
| Rule-based | 0.767 [0.658, 0.868] | 0.958 [0.918, 0.992] | 0.992 [0.974, 1.000] |
| Lookups combined | 0.842 [0.754, 0.924] | 0.958 [0.918, 0.992] | 0.992 [0.974, 1.000] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +0.075 [0.033, 0.126] | +0.000 [0.000, 0.000] | +0.000 [0.000, 0.000] |

#### bg1pct (grpmax below 1%, inheritance filter)

120 cases, median list 14. Rule tie share: Rule-based 0.092.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 0.425 [0.336, 0.517] | 0.658 [0.577, 0.746] | 0.758 [0.667, 0.836] |
| Rule-based | 0.708 [0.598, 0.800] | 0.942 [0.898, 0.976] | 0.942 [0.898, 0.976] |
| Lookups combined | 0.750 [0.655, 0.836] | 0.942 [0.898, 0.976] | 0.942 [0.898, 0.976] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +0.042 [-0.017, 0.104] | +0.000 [0.000, 0.000] | +0.000 [0.000, 0.000] |

#### nomoi (grpmax below 0.1%, no inheritance filter)

120 cases, median list 7. Rule tie share: Rule-based 0.050.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 0.533 [0.449, 0.620] | 0.792 [0.699, 0.874] | 0.992 [0.974, 1.000] |
| Rule-based | 0.742 [0.636, 0.846] | 0.950 [0.910, 0.984] | 0.992 [0.974, 1.000] |
| Lookups combined | 0.817 [0.730, 0.899] | 0.950 [0.910, 0.984] | 0.992 [0.974, 1.000] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +0.075 [0.025, 0.129] | +0.000 [0.000, 0.000] | +0.000 [0.000, 0.000] |

#### all1kg (all of 1000 Genomes left out of the background frequency, inheritance filter)

120 cases, median list 4. Rule tie share: Rule-based 0.042.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| SpliceAI | 0.592 [0.508, 0.686] | 0.892 [0.824, 0.957] | 0.992 [0.974, 1.000] |
| Rule-based | 0.767 [0.661, 0.869] | 0.958 [0.922, 0.992] | 0.992 [0.974, 1.000] |
| Lookups combined | 0.833 [0.754, 0.908] | 0.958 [0.922, 0.992] | 0.992 [0.974, 1.000] |

| Paired contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups combined minus Rule-based | descriptive | +0.067 [0.025, 0.115] | +0.000 [0.000, 0.000] | +0.000 [0.000, 0.000] |

## Per group, top 1 (main analysis)

### SNV

| Score | canonical splice (40) | missense (60) | non-coding (70) | non-coding: intronic (38) | non-coding: splice region (23) | non-coding: synonymous (6) | non-coding: utr (3) | nonsense (30) |
|---|---|---|---|---|---|---|---|---|
| CADD | 0.78 | 0.62 | 0.07 | 0.11 | 0.00 | 0.17 | 0.00 | 0.97 |
| SpliceAI | 0.90 | 0.15 | 0.81 | 0.87 | 0.78 | 1.00 | 0.00 | 0.17 |
| REVEL | 0.03 | 0.90 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| AlphaMissense | 0.00 | 0.78 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| AVI | 0.88 | 0.78 | 0.67 | 0.68 | 0.65 | 0.83 | 0.33 | 0.97 |
| GPN-MSA | 0.55 | 0.60 | 0.30 | 0.37 | 0.30 | 0.00 | 0.00 | 0.40 |
| GPN-Star-V | 0.78 | 0.55 | 0.24 | 0.29 | 0.26 | 0.00 | 0.00 | 0.80 |
| GPN-Star-M | 0.78 | 0.58 | 0.17 | 0.24 | 0.09 | 0.00 | 0.33 | 0.77 |
| Rule-based | 0.97 | 0.85 | 0.60 | 0.63 | 0.61 | 0.67 | 0.00 | 1.00 |
| Rule + AVI | 0.97 | 0.87 | 0.71 | 0.79 | 0.70 | 0.67 | 0.00 | 1.00 |
| Lookups combined | 0.82 | 0.85 | 0.77 | 0.89 | 0.70 | 0.67 | 0.00 | 1.00 |
| Lookups + AVI | 0.85 | 0.83 | 0.89 | 0.97 | 0.87 | 0.83 | 0.00 | 1.00 |
| Lookups + GPN-MSA | 0.85 | 0.83 | 0.73 | 0.87 | 0.65 | 0.50 | 0.00 | 0.97 |
| Lookups + GPN-Star-V | 0.82 | 0.85 | 0.77 | 0.92 | 0.65 | 0.67 | 0.00 | 0.97 |
| Lookups + GPN-Star-M | 0.85 | 0.85 | 0.79 | 0.89 | 0.78 | 0.50 | 0.00 | 0.97 |

Groups with fewer than 20 cases are small: synonymous (6) and UTR (3).

Primary and key secondary by group (top 1 gain, 95% paired interval):

| Group | Lookups + AVI minus Lookups combined | Rule + AVI minus Rule-based |
|---|---|---|
| canonical splice | +0.025 [-0.125, 0.171] | +0.000 [-0.071, 0.075] |
| missense | -0.017 [-0.108, 0.066] | +0.017 [0.000, 0.053] |
| non-coding | +0.114 [0.027, 0.224] | +0.114 [0.025, 0.216] |
| non-coding: intronic | +0.079 [-0.029, 0.220] | +0.158 [0.049, 0.297] |
| non-coding: splice region | +0.174 [0.040, 0.364] | +0.087 [0.000, 0.200] |
| non-coding: synonymous | +0.167 [0.000, 0.500] | +0.000 [-0.500, 0.500] |
| non-coding: utr | +0.000 [0.000, 0.000] | +0.000 [0.000, 0.000] |
| nonsense | +0.000 [0.000, 0.000] | +0.000 [0.000, 0.000] |

Both AVI gains come from the non-coding group (+0.114 each). Lookups + AVI gains most on splice region diagnostics, Rule + AVI on intronic ones. Neither gains on nonsense; the gains on canonical splice and missense are near 0 with wide intervals.

### Indel

| Score | canonical splice (30) | frameshift (40) | inframe (30) | non-coding (20) | non-coding: intronic (13) | non-coding: utr (7) |
|---|---|---|---|---|---|---|
| SpliceAI | 1.00 | 0.33 | 0.47 | 0.65 | 0.92 | 0.14 |
| Rule-based | 1.00 | 0.88 | 0.70 | 0.30 | 0.38 | 0.14 |
| Lookups combined | 1.00 | 0.93 | 0.77 | 0.55 | 0.77 | 0.14 |

Groups with fewer than 20 cases are small: intronic (13) and UTR (7). Six of the seven UTR diagnostics have no SpliceAI score, so every score ranks them last (0.14). The rule gives intronic indels only a SpliceAI level, so it ranks them below Lookups combined (0.38 against 0.77).

## Diagnostic found in gnomAD or absent (section 9)

The diagnostic's status comes from `found` in data/gnomad_full.tsv, then data/pool/gnomad_full.tsv. Every diagnostic has a status. SNV: 53 found, 147 absent. Indel: 29 found, 91 absent. Both contrasts in each stratum are on one shared draw of people (seed "<study>/<analysis>/gnomad split"), so the difference absent minus found is paired.

| SNV analysis | Contrast | Found | Absent | Absent minus found |
|---|---|---|---|---|
| main | primary | +0.057 [-0.020, 0.149] | +0.034 [-0.028, 0.101] | -0.023 [-0.119, 0.060] |
| main | key secondary | +0.000 [-0.065, 0.050] | +0.061 [0.020, 0.116] | +0.061 [-0.006, 0.138] |
| bg1pct | primary | +0.038 [-0.127, 0.184] | +0.000 [-0.094, 0.094] | -0.038 [-0.189, 0.132] |
| bg1pct | key secondary | +0.000 [-0.046, 0.052] | +0.054 [0.007, 0.105] | +0.054 [-0.011, 0.127] |
| nomoi | primary | +0.094 [0.000, 0.200] | +0.041 [-0.047, 0.134] | -0.053 [-0.167, 0.064] |
| nomoi | key secondary | +0.038 [-0.038, 0.115] | +0.082 [0.029, 0.138] | +0.044 [-0.038, 0.133] |
| all1kg | primary | +0.075 [-0.020, 0.184] | +0.075 [0.007, 0.159] | -0.001 [-0.120, 0.121] |
| all1kg | key secondary | -0.019 [-0.067, 0.000] | +0.068 [0.026, 0.124] | +0.087 [0.033, 0.156] |

- **Primary:** `possibly_gnomad_presence` is false in every analysis. The gain from AVI in Lookups + AVI is no larger in absent cases; in main it is slightly larger in found cases.
- **Key secondary:** the code flags the primary only, but section 9 asks for a gain only in absent cases to be reported. Rule + AVI shows that pattern. Its gain is 0 or below in found cases in main, bg1pct and all1kg, and positive with an interval above 0 in absent cases in every analysis. The difference excludes 0 in all1kg (+0.087 [0.033, 0.156]) and nearly does in main (+0.061 [-0.006, 0.138]). So the key secondary gain is possibly due to gnomAD presence, not to what a lab would see.
- Top 1 by stratum, main: AVI alone 0.717 in found cases against 0.816 in absent cases; CADD 0.509 against 0.510; Rule-based 0.792 against 0.816. AVI ranks absent diagnostics higher than found ones, while CADD does not.
- Indel: no AVI contrast, so the split has cases only. Top 1 in main, found against absent: SpliceAI 0.655 and 0.560, Rule-based 0.828 and 0.747, Lookups combined 0.862 and 0.835.

## AVI training chromosomes (snv heldout)

Cases split by whether the diagnostic is on an AVI training chromosome (1, 4, 7, 8, 10, 13, 15); each ranked only among background variants in the same split.

| Split | Cases | Median list | AVI | Lookups combined | Lookups + AVI | Rule-based | Rule + AVI |
|---|---|---|---|---|---|---|---|
| training | 38 | 7 | 0.816 | 0.947 | 0.974 | 0.895 | 0.921 |
| other | 162 | 19 | 0.840 | 0.883 | 0.914 | 0.827 | 0.870 |

- Training chromosomes: primary +0.026 [-0.063, 0.108]; key secondary +0.026 [0.000, 0.089].
- Other chromosomes: primary +0.031 [-0.006, 0.075]; key secondary +0.043 [0.006, 0.082].
- Difference, training minus other: primary -0.004 [-0.098, 0.086]; key secondary -0.017 [-0.075, 0.056]. `possibly_inflated_by_training_overlap` is false.
- Previously the primary difference was -0.102 [-0.190, -0.024]. The change is the Lookups combined refit described above, on a split of 38 cases with a median list of 7: Lookups combined went from 0.974 to 0.947 there and Lookups + AVI from 0.921 to 0.974. One or two cases move these numbers by 3 to 5 points.

## EVEE (snv evee, descriptive only)

Agreement with ClinVar labels on in-sample variants, not diagnostic performance (section 10). bg1pct filters, lists restricted to variants with an EVEE row (7,084 of the 19,561). Every one of the 200 diagnostics has a row, so no case was dropped; median list 29.5. Not to be compared with the main analysis.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 0.445 [0.375, 0.515] | 0.740 [0.680, 0.805] | 0.910 [0.860, 0.955] |
| SpliceAI | 0.500 [0.435, 0.570] | 0.635 [0.580, 0.685] | 0.660 [0.610, 0.710] |
| REVEL | 0.265 [0.215, 0.315] | 0.315 [0.260, 0.370] | 0.330 [0.275, 0.385] |
| AlphaMissense | 0.240 [0.195, 0.285] | 0.280 [0.235, 0.325] | 0.285 [0.240, 0.330] |
| AVI | 0.745 [0.665, 0.815] | 0.935 [0.895, 0.965] | 0.975 [0.950, 0.995] |
| GPN-MSA | 0.415 [0.335, 0.490] | 0.745 [0.695, 0.790] | 0.870 [0.830, 0.910] |
| GPN-Star-V | 0.515 [0.440, 0.585] | 0.820 [0.775, 0.870] | 0.920 [0.880, 0.955] |
| GPN-Star-M | 0.425 [0.340, 0.510] | 0.810 [0.760, 0.855] | 0.920 [0.885, 0.955] |
| EVEE | 0.815 [0.755, 0.885] | 0.975 [0.950, 0.995] | 0.985 [0.965, 1.000] |
| Lookups combined | 0.820 [0.745, 0.885] | 0.965 [0.930, 0.990] | 0.975 [0.940, 0.995] |

EVEE and Lookups combined are close (0.815 and 0.820), as before (0.815 and 0.825). The small changes come from the larger bg1pct lists (19,309 to 19,561 variants; 7,004 to 7,084 with an EVEE row) and the refit.

## Checks

- **Invariance.** `tests/test_invariance.py` with GPN-Star-M dropped, results written to a scratch folder: snv main 1,707 rows compared, 267 involve GPN-Star-M or are skipped, 0 problems; indel main 348 rows compared, 1 skipped, 0 problems.
- **Lookups combined against single tools.** It beats every single tool overall in main (0.840 against AVI 0.790), bg1pct (0.770 against 0.700), nomoi (0.685 against 0.665) and evee (0.820 against EVEE 0.815), but not in all1kg (0.770 against AVI 0.775). It beats SpliceAI in every indel analysis. Within groups it falls below a single tool in places: canonical splice 0.82 against SpliceAI 0.90 and AVI 0.88, and indel intronic 0.77 against SpliceAI 0.92 (13 cases).
- **Charts.** `figures/rd_charts.py` wrote gains_primary, gains_added, gains_gnomad and top1_by_group for snv and indel to `figures/out/`, from these files. There are no replication results yet, so the charts hold the cross-validated rows only. One layout fault: panels in the top row with no panel beneath them (indel canonical splice, frameshift, inframe and non-coding; SNV non-coding: intronic) show no x tick numbers, because the shared x axis keeps tick labels on the bottom row only.
