# Results on the 500 replication cases (locked)

These are the primary results of the second run (PROTOCOL.md section 8). The 500 cases come from 50 people that the first run never used, with 10 cases each. Every fitted model was fitted once on the 200 SNV cases under the main filters and applied unchanged here ("locked"); the 1% filter check and EVEE use models fitted on the 200 under the 1% filter. Rule-based and Rule + AVI have nothing to fit.

Made by `case_benchmark.py --study replication` for main, bg1pct, nomoi, all1kg, heldout and evee: `results/replication_<analysis>.json`, logs in `logs/bench_replication_<analysis>.log`. Package versions: python 3.12.14, numpy 2.5.4, pandas 3.1.0, scikit-learn 1.9.1, scipy 1.18.1, joblib 1.6.0, threadpoolctl 3.7.0.

All figures are percentages, rounded to one decimal. Gains are in percentage points. Intervals are 95% bootstrap intervals over people (1,000 resamples, shared by every score in a stratum), and every gain has a paired interval; the GPN family also has a Bonferroni-adjusted 98.33% interval. Top k is the share of cases with fewer than k competitors: background variants that score at least as high as the diagnostic (for the rule scores, background variants whose pairwise comparison with the diagnostic is TIED or BEATS, section 7). Top 1 is the share of cases where the diagnostic ranks first.

Main analysis: background grpmax AF below 0.1% with the person and their family left out of gnomAD, and the background inheritance filter. 500 cases, 50 people, median list 19 SNVs (diagnostic included).

## Primary result

Adding AVI to the lookups combined raises top 1 by +1.8 points (-0.4 to +4.0), from 90.6% to 92.4%. The interval lies inside -5 to +5, so by the reading rule fixed in advance the verdict is **no gain that matters**. `case_benchmark.py` writes the same verdict (`primary.verdict`).

Key secondary: Rule + AVI minus Rule-based is +2.6 points (+0.4 to +5.6), from 87.2% to 89.8%. The interval is above 0 but its upper end is past +5, so read by the same rule on its own it would be inconclusive. It does not change the primary verdict: the rule uses the key secondary only for "a gain that matters", which the primary does not reach.

| Contrast (top 1) | Role | Gain | Reading |
|---|---|---|---|
| Lookups + AVI minus Lookups combined, locked | primary | +1.8 (-0.4 to +4.0) | No gain that matters |
| Rule + AVI minus Rule-based | key secondary | +2.6 (+0.4 to +5.6) | above 0, estimate under 5 points, interval reaches +5.6 |
| Lookups + AVI minus Lookups combined, cross-validated on the 500 | descriptive | +0.8 (-1.8 to +3.2) | |

At top 5 the primary gain is +1.0 (+0.2 to +2.2) and at top 10 +0.2 (-0.4 to +0.8).

## All scores, all cases

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 61.0 (53.4 to 68.8) | 91.2 (87.6 to 94.4) | 97.6 (95.4 to 99.4) |
| SpliceAI | 57.8 (53.0 to 62.6) | 70.4 (66.8 to 74.4) | 79.4 (75.0 to 83.8) |
| REVEL | 28.6 (24.4 to 32.4) | 36.8 (31.8 to 42.4) | 47.6 (40.0 to 56.2) |
| AlphaMissense | 24.2 (20.6 to 28.0) | 33.0 (27.8 to 39.0) | 44.0 (35.8 to 52.8) |
| AVI | 83.6 (78.2 to 88.8) | 98.0 (96.8 to 99.0) | 99.6 (99.0 to 100.0) |
| GPN-MSA | 51.2 (43.8 to 58.2) | 83.8 (79.4 to 88.0) | 93.2 (90.4 to 95.8) |
| GPN-Star-V | 59.4 (51.8 to 66.0) | 87.0 (82.6 to 90.8) | 95.6 (93.4 to 97.4) |
| GPN-Star-M | 55.8 (48.0 to 63.4) | 87.2 (83.4 to 90.8) | 95.0 (92.8 to 97.0) |
| Rule-based | 87.2 (82.8 to 91.0) | 97.8 (96.6 to 98.8) | 98.6 (97.6 to 99.6) |
| Rule + AVI | 89.8 (86.6 to 92.8) | 98.0 (96.8 to 99.2) | 98.6 (97.6 to 99.4) |
| Lookups combined (locked) | 90.6 (87.8 to 93.4) | 97.0 (95.4 to 98.4) | 98.6 (97.4 to 99.6) |
| Lookups + AVI (locked) | 92.4 (89.8 to 94.8) | 98.0 (96.6 to 99.2) | 98.8 (97.4 to 99.8) |
| Lookups + GPN-MSA (locked) | 90.4 (87.4 to 93.2) | 96.8 (95.2 to 98.2) | 98.4 (97.2 to 99.4) |
| Lookups + GPN-Star-V (locked) | 91.0 (87.8 to 93.6) | 98.0 (96.8 to 99.2) | 98.4 (97.2 to 99.4) |
| Lookups + GPN-Star-M (locked) | 90.4 (87.4 to 93.4) | 97.0 (95.2 to 98.4) | 98.2 (96.8 to 99.2) |
| Lookups combined (cross-validated on the 500) | 92.8 (89.8 to 95.6) | 98.4 (97.4 to 99.4) | 99.4 (98.6 to 100.0) |
| Lookups + AVI (cross-validated on the 500) | 93.6 (90.6 to 96.2) | 98.6 (97.6 to 99.6) | 99.8 (99.4 to 100.0) |

The lookups are near the ceiling. Locked, they rank the diagnostic first in 90.6% (87.8 to 93.4) of cases and in the top 5 in 97.0% (95.4 to 98.4). AVI alone ranks it first in 83.6% (78.2 to 88.8), more than any other single score (next: CADD 61.0%, 53.4 to 68.8; GPN-Star-V 59.4%, 51.8 to 66.0), but less than either baseline.

## Secondary: the GPN family (Bonferroni-adjusted)

Each GPN score added to the lookups combined, locked, top 1. A secondary score is said to add something only if its 98.33% interval excludes zero (section 8).

| Added score | Top 1 gain | 95% interval | 98.33% interval | Adds something? |
|---|---|---|---|---|
| GPN-MSA | -0.2 | (-1.4 to +1.0) | (-1.6 to +1.4) | No |
| GPN-Star-V | +0.4 | (-1.2 to +2.0) | (-1.8 to +2.4) | No |
| GPN-Star-M | -0.2 | (-1.4 to +1.0) | (-1.6 to +1.2) | No |

None of the three adds anything. In the first run, GPN-MSA's +6 point gain on the 200 cases also failed to replicate on its 500.

## Descriptive contrasts, all cases

Unadjusted 95% intervals, descriptive (section 8).

| Contrast | Role | Top 1 gain | Top 5 gain | Top 10 gain |
|---|---|---|---|---|
| Lookups + AVI minus Lookups combined (locked) | primary | +1.8 (-0.4 to +4.0) | +1.0 (+0.2 to +2.2) | +0.2 (-0.4 to +0.8) |
| Rule + AVI minus Rule-based | key secondary | +2.6 (+0.4 to +5.6) | +0.2 (0.0 to +0.6) | 0.0 (-0.6 to +0.6) |
| Lookups combined (locked) minus Rule-based | descriptive | +3.4 (0.0 to +7.4) | -0.8 (-2.4 to +0.8) | 0.0 (-1.2 to +1.2) |
| Lookups + AVI (locked) minus Rule-based | descriptive | +5.2 (+1.6 to +9.2) | +0.2 (-1.0 to +1.6) | +0.2 (-1.2 to +1.4) |
| Lookups + AVI minus Lookups combined (cross-validated) | primary contrast, cross-validated (descriptive) | +0.8 (-1.8 to +3.2) | +0.2 (0.0 to +0.6) | +0.4 (0.0 to +1.0) |
| Lookups combined (cross-validated) minus Rule-based | descriptive | +5.6 (+1.2 to +10.0) | +0.6 (-0.4 to +1.6) | +0.8 (-0.2 to +2.0) |
| Lookups + AVI (cross-validated) minus Rule-based | descriptive | +6.4 (+1.8 to +11.0) | +0.8 (-0.2 to +2.0) | +1.2 (+0.2 to +2.2) |

The locked lookups combined beat the rule-based ranking by +3.4 (0.0 to +7.4) points at top 1, on the edge of zero. By group (below), the lookups lead on non-coding diagnostics (+12.0, +4.7 to +20.9) and trail on missense (-2.7, -7.3 to +2.8). At top 5 the overall difference is -0.8 (-2.4 to +0.8).

## By diagnostic group

Top 1 by group, and by CADD class within the non-coding group. Lookups combined and Lookups + AVI are the locked models. Groups under 20 cases are marked small; their intervals are wide.

| Group | Cases | Lookups combined (locked) | Lookups + AVI (locked) | Rule-based | Rule + AVI | AVI |
|---|---|---|---|---|---|---|
| Canonical splice | 100 | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 98.0 (94.7 to 100.0) |
| Missense | 150 | 87.3 (82.1 to 92.5) | 87.3 (80.8 to 93.0) | 90.0 (84.0 to 95.4) | 88.7 (82.8 to 93.9) | 80.7 (72.0 to 89.4) |
| Non-coding | 175 | 84.0 (77.4 to 89.9) | 89.1 (85.0 to 93.3) | 72.0 (61.7 to 81.1) | 80.6 (72.7 to 87.6) | 71.4 (62.2 to 80.2) |
| Non-coding: intronic | 84 | 89.3 (81.2 to 96.3) | 94.0 (89.2 to 98.7) | 76.2 (62.3 to 87.5) | 89.3 (80.8 to 95.7) | 72.6 (60.7 to 83.3) |
| Non-coding: splice region | 63 | 93.7 (87.7 to 98.4) | 92.1 (85.0 to 97.4) | 77.8 (63.1 to 89.2) | 81.0 (66.7 to 91.9) | 77.8 (63.8 to 88.9) |
| Non-coding: synonymous | 15 (small) | 86.7 (66.7 to 100.0) | 100.0 (100.0 to 100.0) | 86.7 (66.7 to 100.0) | 86.7 (66.7 to 100.0) | 73.3 (50.0 to 93.8) |
| Non-coding: UTR | 13 (small) | 0.0 (0.0 to 0.0) | 30.8 (8.3 to 50.0) | 0.0 (0.0 to 0.0) | 15.4 (0.0 to 33.3) | 30.8 (8.3 to 50.0) |
| Nonsense | 75 | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 98.7 (95.9 to 100.0) |

Gains by group, top 1, 95% paired interval (descriptive):

| Group | Cases | Primary (Lookups + AVI minus Lookups combined, locked) | Key secondary (Rule + AVI minus Rule-based) | Lookups combined (locked) minus Rule-based |
|---|---|---|---|---|
| Canonical splice | 100 | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) |
| Missense | 150 | 0.0 (-5.5 to +4.6) | -1.3 (-3.8 to +1.3) | -2.7 (-7.3 to +2.8) |
| Non-coding | 175 | +5.1 (+0.6 to +10.5) | +8.6 (+3.2 to +15.1) | +12.0 (+4.7 to +20.9) |
| Non-coding: intronic | 84 | +4.8 (0.0 to +11.6) | +13.1 (+3.2 to +25.0) | +13.1 (+2.1 to +26.8) |
| Non-coding: splice region | 63 | -1.6 (-8.8 to +5.4) | +3.2 (0.0 to +8.8) | +15.9 (+3.2 to +28.6) |
| Non-coding: synonymous | 15 | +13.3 (0.0 to +33.3) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) |
| Non-coding: UTR | 13 | +30.8 (+8.3 to +50.0) | +15.4 (0.0 to +33.3) | 0.0 (0.0 to 0.0) |
| Nonsense | 75 | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) |

- Canonical splice and nonsense diagnostics rank first in every case under both baselines, with or without AVI. Nothing can be gained there.
- Missense: the primary gain is 0.0 (-5.5 to +4.6), and the rule-based ranking is as good as the lookups (90.0%, 84.0 to 95.4, against 87.3%, 82.1 to 92.5).
- Non-coding is where AVI adds: primary +5.1 (+0.6 to +10.5), key secondary +8.6 (+3.2 to +15.1). The rule-based ranking is weakest here (72.0%, 61.7 to 81.1, against 84.0%, 77.4 to 89.9, for the lookups), because outside missense it has only SpliceAI's supporting level to go on.
- UTR diagnostics (13 cases) are never ranked first by either baseline; AVI alone and Lookups + AVI rank 30.8% (8.3 to 50.0) first. Synonymous (15) and UTR are too small to read on their own.

Single scores by group, top 1 (descriptive):

| Group | Cases | CADD | SpliceAI | REVEL | AlphaMissense | AVI | GPN-MSA | GPN-Star-V | GPN-Star-M |
|---|---|---|---|---|---|---|---|---|---|
| Canonical splice | 100 | 89.0 (80.4 to 95.6) | 96.0 (91.1 to 100.0) | 1.0 (0.0 to 3.2) | 0.0 (0.0 to 0.0) | 98.0 (94.7 to 100.0) | 72.0 (57.8 to 84.0) | 81.0 (70.6 to 90.3) | 79.0 (70.5 to 87.9) |
| Missense | 150 | 63.3 (50.7 to 75.2) | 14.7 (9.1 to 20.9) | 90.7 (84.7 to 95.7) | 80.7 (72.5 to 87.5) | 80.7 (72.0 to 89.4) | 64.0 (53.5 to 74.0) | 69.3 (59.6 to 79.5) | 60.0 (48.0 to 71.9) |
| Non-coding | 175 | 26.9 (16.3 to 38.2) | 86.3 (81.3 to 91.3) | 1.7 (0.0 to 3.9) | 0.0 (0.0 to 0.0) | 71.4 (62.2 to 80.2) | 28.6 (19.9 to 37.7) | 32.0 (22.6 to 41.7) | 30.9 (20.8 to 41.5) |
| Non-coding: intronic | 84 | 27.4 (16.0 to 40.2) | 94.0 (88.9 to 98.8) | 2.4 (0.0 to 6.1) | 0.0 (0.0 to 0.0) | 72.6 (60.7 to 83.3) | 34.5 (22.1 to 46.1) | 29.8 (19.3 to 40.0) | 28.6 (18.9 to 37.8) |
| Non-coding: splice region | 63 | 25.4 (10.8 to 40.6) | 90.5 (81.8 to 96.8) | 1.6 (0.0 to 5.3) | 0.0 (0.0 to 0.0) | 77.8 (63.8 to 88.9) | 22.2 (8.9 to 35.7) | 33.3 (17.6 to 48.5) | 30.2 (14.3 to 46.0) |
| Non-coding: synonymous | 15 (small) | 33.3 (13.3 to 56.3) | 93.3 (78.6 to 100.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 73.3 (50.0 to 93.8) | 20.0 (0.0 to 42.9) | 20.0 (0.0 to 42.9) | 33.3 (11.1 to 57.1) |
| Non-coding: UTR | 13 (small) | 23.1 (0.0 to 46.2) | 7.7 (0.0 to 25.0) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) | 30.8 (8.3 to 50.0) | 30.8 (7.7 to 58.3) | 53.8 (30.8 to 76.9) | 46.2 (23.1 to 69.2) |
| Nonsense | 75 | 98.7 (95.5 to 100.0) | 26.7 (17.1 to 36.4) | 4.0 (0.0 to 9.2) | 0.0 (0.0 to 0.0) | 98.7 (95.9 to 100.0) | 50.7 (38.0 to 63.2) | 74.7 (63.1 to 84.7) | 74.7 (61.1 to 86.4) |

## Checks (section 9)

### Summary

The reading rule applies to the main analysis only. `case_benchmark.py` writes a verdict string for every analysis; for the checks it is shown for reference and is not a result.

| Analysis | Median list | Lookups combined (locked) | Lookups + AVI (locked) | Rule-based | Rule + AVI | AVI | Primary gain | Key secondary gain | Verdict string |
|---|---|---|---|---|---|---|---|---|---|
| main (grpmax below 0.1%, family left out, inheritance filter) | 19 | 90.6 (87.8 to 93.4) | 92.4 (89.8 to 94.8) | 87.2 (82.8 to 91.0) | 89.8 (86.6 to 92.8) | 83.6 (78.2 to 88.8) | +1.8 (-0.4 to +4.0) | +2.6 (+0.4 to +5.6) | No gain that matters |
| bg1pct (grpmax below 1%; locked models fitted under 1%) | 69.5 | 82.0 (76.4 to 87.2) | 80.8 (74.0 to 86.4) | 74.8 (68.0 to 81.0) | 80.0 (73.8 to 85.4) | 71.4 (65.2 to 77.4) | -1.2 (-6.6 to +3.2) | +5.2 (+1.4 to +9.6) | Inconclusive |
| nomoi (no background inheritance filter) | 34.5 | 79.6 (73.2 to 85.4) | 83.6 (77.4 to 89.2) | 72.4 (65.2 to 79.2) | 74.6 (67.6 to 80.8) | 71.2 (63.6 to 78.2) | +4.0 (0.0 to +8.0) | +2.2 (-1.2 to +5.4) | Inconclusive |
| all1kg (all of 1000 Genomes left out of background frequencies) | 21 | 87.8 (83.2 to 91.8) | 88.2 (83.0 to 92.6) | 84.2 (78.6 to 89.2) | 86.8 (81.6 to 91.2) | 80.4 (74.2 to 85.6) | +0.4 (-2.8 to +3.0) | +2.6 (+0.4 to +5.4) | No gain that matters |

- **Looser background filter (bg1pct).** Lists grow to a median of 69.5 and every score falls. AVI's gain to the lookups is -1.2 (-6.6 to +3.2); its gain to the rule-based ranking is +5.2 (+1.4 to +9.6).
- **No background inheritance filter (nomoi).** Lists grow to a median of 34.5. The primary gain is +4.0 (0.0 to +8.0), the largest of any analysis, and still under 5 points; the key secondary is +2.2 (-1.2 to +5.4).
- **All of 1000 Genomes left out (all1kg).** Lists barely change (median 21). Primary +0.4 (-2.8 to +3.0), key secondary +2.6 (+0.4 to +5.4). Leaving the family out and leaving all of 1000 Genomes out give the same answer.

GPN family under each check, top 1 gain with the 98.33% interval:

| Analysis | GPN-MSA | GPN-Star-V | GPN-Star-M |
|---|---|---|---|
| main | -0.2 (-1.6 to +1.4) | +0.4 (-1.8 to +2.4) | -0.2 (-1.6 to +1.2) |
| bg1pct | +0.4 (-2.4 to +3.6) | +0.4 (-3.0 to +3.8) | +2.2 (-0.2 to +5.1) |
| nomoi | -0.2 (-3.2 to +2.2) | -0.2 (-3.7 to +3.2) | +1.4 (-1.6 to +5.3) |
| all1kg | -1.4 (-3.4 to 0.0) | +0.6 (-1.6 to +2.8) | 0.0 (-1.4 to +1.4) |

No 98.33% interval excludes zero. The nearest are GPN-Star-M under bg1pct, whose interval starts just below zero, and GPN-MSA under all1kg, whose interval ends at zero (a loss, if anything).

### AVI training chromosomes (heldout)

Cases split by whether the diagnostic is on an AVI training chromosome (1, 4, 7, 8, 10, 13, 15); each case is ranked only among its background variants on chromosomes of the same split, with the locked models unchanged. Both splits use the same draws of people, so the difference is paired.

| Split | Cases | Median list | AVI | Lookups combined (locked) | Lookups + AVI (locked) | Rule-based | Rule + AVI | Primary gain | Key secondary gain |
|---|---|---|---|---|---|---|---|---|---|
| training chromosomes | 124 | 5 | 90.3 (82.0 to 97.0) | 93.5 (88.1 to 98.2) | 95.2 (90.2 to 99.2) | 95.2 (91.4 to 98.3) | 96.0 (92.2 to 99.2) | +1.6 (0.0 to +4.2) | +0.8 (-2.6 to +4.6) |
| other chromosomes | 376 | 13 | 86.2 (80.6 to 91.2) | 91.2 (88.1 to 94.1) | 93.3 (90.6 to 96.0) | 88.3 (84.3 to 92.0) | 90.7 (86.8 to 93.9) | +2.1 (-0.3 to +4.6) | +2.4 (+0.5 to +4.7) |

Difference, training minus other: primary -0.5 (-3.7 to +3.1), key secondary -1.6 (-5.6 to +2.5). The gain is not larger on AVI's training chromosomes, so the primary result is not flagged (`possibly_inflated_by_training_overlap` is false).

### Diagnostic found in gnomAD or absent

Found means the diagnostic is in full gnomAD v3.1.2 genomes with AC above 0: 114 found, 386 absent, 0 without a status. Both contrasts in each analysis use one shared draw of people, so absent minus found is paired. A gain is flagged (`possibly_gnomad_presence`) if it is larger in the absent cases and that interval excludes zero.

| Analysis | Contrast | Found | Absent | Absent minus found | Flagged |
|---|---|---|---|---|---|
| main | primary | -0.9 (-6.6 to +5.0) | +2.6 (+0.5 to +4.8) | +3.5 (-2.2 to +8.8) | No |
| main | key secondary | +3.5 (0.0 to +8.2) | +2.3 (0.0 to +5.2) | -1.2 (-5.3 to +2.8) | No |
| bg1pct | primary | -2.6 (-9.6 to +4.1) | -0.8 (-6.2 to +4.2) | +1.8 (-4.1 to +7.8) | No |
| bg1pct | key secondary | +8.8 (+2.6 to +16.2) | +4.2 (0.0 to +8.4) | -4.6 (-10.1 to +0.8) | No |
| nomoi | primary | +3.5 (-3.5 to +11.0) | +4.2 (+0.3 to +8.6) | +0.6 (-6.4 to +7.2) | No |
| nomoi | key secondary | +0.9 (-6.5 to +7.0) | +2.6 (-0.3 to +6.2) | +1.7 (-4.6 to +9.0) | No |
| all1kg | primary | -2.6 (-9.7 to +4.4) | +1.3 (-1.0 to +3.5) | +3.9 (-2.4 to +9.8) | No |
| all1kg | key secondary | +4.4 (+0.8 to +9.4) | +2.1 (-0.3 to +4.9) | -2.3 (-6.8 to +2.1) | No |

No gain is flagged in any analysis. On the 500 the key secondary is no larger in absent cases than in found ones, unlike the 200 (see `analysis_snv_indel.md`).

Top 1 in main, found against absent: AVI 75.4 (65.6 to 84.2) against 86.0 (80.0 to 91.2); CADD 56.1 (45.4 to 67.0) against 62.4 (54.2 to 70.2); Lookups combined (locked) 89.5 (84.0 to 94.7) against 90.9 (87.7 to 93.9); Rule-based 83.3 (76.3 to 89.8) against 88.3 (83.4 to 93.0).

## EVEE (descriptive only)

This is agreement with ClinVar labels on in-sample variants, not diagnostic performance, and it is not to be compared with the main analysis (section 10). EVEE's scores cover ClinVar variants only, almost every diagnostic was already P/LP in EVEE's training labels, and most covered background variants are B/LB, its training negatives. Whether the released scores are held-out predictions is not documented.

bg1pct filters (models fitted on the 200 under the 1% filter, on the full lists), then each list restricted to variants with an EVEE row. 8 of the 500 diagnostics have no EVEE row, so 492 cases remain; median list 26; 15,398 variants in the lists.

| Score | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| CADD | 51.0 (44.5 to 58.1) | 83.7 (78.8 to 88.3) | 92.5 (89.4 to 95.4) |
| SpliceAI | 51.6 (46.2 to 56.4) | 64.0 (59.7 to 68.0) | 67.5 (62.7 to 71.5) |
| REVEL | 26.2 (22.5 to 29.8) | 31.7 (28.0 to 35.5) | 32.7 (28.8 to 36.8) |
| AlphaMissense | 24.2 (20.8 to 27.5) | 28.7 (24.9 to 32.1) | 29.1 (25.4 to 32.6) |
| AVI | 76.0 (69.7 to 82.0) | 95.5 (93.5 to 97.4) | 98.8 (97.8 to 99.6) |
| GPN-MSA | 44.5 (37.1 to 51.7) | 78.9 (74.8 to 82.5) | 90.5 (87.4 to 93.3) |
| GPN-Star-V | 52.2 (45.7 to 58.4) | 83.3 (78.5 to 87.6) | 92.9 (90.4 to 95.3) |
| GPN-Star-M | 49.4 (43.5 to 55.5) | 80.1 (75.4 to 84.6) | 94.5 (92.5 to 96.4) |
| EVEE | 87.0 (82.5 to 90.9) | 99.2 (98.4 to 99.8) | 99.8 (99.4 to 100.0) |
| Lookups combined (locked) | 87.0 (82.2 to 91.0) | 95.7 (93.7 to 97.6) | 96.8 (94.7 to 98.6) |
| Lookups combined (cross-validated on the 500) | 88.0 (82.8 to 92.0) | 97.8 (96.5 to 98.8) | 98.6 (97.5 to 99.4) |

EVEE alone matches the locked lookups combined at top 1 (87.0%, 82.5 to 90.9, and 87.0%, 82.2 to 91.0). Given its training labels, that is what agreement with ClinVar would look like; it says nothing about how EVEE would rank a new patient's variants.

| Group | Cases | EVEE | Lookups combined (locked) | AVI |
|---|---|---|---|---|
| Canonical splice | 98 | 99.0 (96.7 to 100.0) | 92.9 (86.4 to 98.0) | 96.9 (92.8 to 100.0) |
| Missense | 150 | 80.7 (71.9 to 88.5) | 84.0 (74.0 to 92.4) | 74.7 (65.1 to 84.7) |
| Non-coding | 172 | 80.8 (74.1 to 86.7) | 81.4 (75.2 to 87.4) | 55.8 (46.8 to 65.4) |
| Non-coding: intronic | 81 | 77.8 (69.6 to 86.6) | 87.6 (79.7 to 94.5) | 59.3 (46.2 to 71.6) |
| Non-coding: splice region | 63 | 90.5 (83.0 to 97.0) | 88.9 (80.4 to 96.7) | 61.9 (48.1 to 75.0) |
| Non-coding: synonymous | 15 (small) | 73.3 (53.3 to 92.9) | 86.7 (64.3 to 100.0) | 60.0 (33.3 to 83.3) |
| Non-coding: UTR | 13 (small) | 61.5 (40.0 to 83.3) | 0.0 (0.0 to 0.0) | 0.0 (0.0 to 0.0) |
| Nonsense | 72 | 98.6 (95.2 to 100.0) | 98.6 (94.9 to 100.0) | 98.6 (95.3 to 100.0) |

## List lengths

Median list length per case, diagnostic included. No case lost its diagnostic under any filter (`cases_lost_diagnostic` is 0 in every run).

| Analysis | Cases | Median list | SNVs ranked | Rows after the frequency filter | Rows after the inheritance filter (both types) |
|---|---|---|---|---|---|
| main | 500 | 19 | 10,475 | 20,580 | 11,856 |
| bg1pct | 500 | 69.5 | 42,321 | 80,870 | 49,512 |
| nomoi | 500 | 34.5 | 18,620 | 20,580 | not applied |
| all1kg | 500 | 21 | 10,985 | 21,670 | 12,486 |
| evee | 492 | 26 | 15,398 | 80,870 | 49,512 |
| heldout | 500 | 5 (training split), 13 (other) | | | |

For evee, SNVs ranked counts only variants with an EVEE row. In main, 1,381 background indels count toward the two-variant rule in biallelic genes and are then dropped, as the study ranks SNVs.

For comparison, the first run's grpmax side check left a median of 34.5 background variants on its 500 cases (`../rare-disease/REPLICATION.md`). Here the same frequency filter without the inheritance filter (nomoi) leaves a median list of 34.5, diagnostic included; the inheritance filter brings it down to 19.

## Rule tie share

Share of cases where at least one background variant TIES the diagnostic under the pairwise rule (section 7).

| Analysis | Rule-based | Rule + AVI |
|---|---|---|
| main | 0.0% (0 of 500) | 0.0% (0 of 500) |
| bg1pct | 0.4% (2 of 500) | 0.0% (0 of 500) |
| nomoi | 0.6% (3 of 500) | 0.0% (0 of 500) |
| all1kg | 0.0% (0 of 500) | 0.0% (0 of 500) |

The few ties are non-coding diagnostics in bg1pct and nomoi, where a background variant gets the same points from the same tool and the rule cannot separate them. Rule + AVI has none, because AVI breaks them.

## Coverage (main)

Variants with a score, over all variants in the lists (diagnostics included), and diagnostics with a score.

| Score | Variants (of 10,475) | Diagnostics (of 500) |
|---|---|---|
| CADD | 10,475 | 500 |
| SpliceAI | 8,729 | 484 |
| REVEL | 2,090 | 168 |
| AlphaMissense | 1,885 | 143 |
| AVI | 10,475 | 500 |
| GPN-MSA | 10,475 | 500 |
| GPN-Star-V | 10,475 | 500 |
| GPN-Star-M | 10,475 | 500 |

REVEL and AlphaMissense score missense variants only. A missing score ranks last.

## Model classes

From the 200 SNV cases: canonical splice 40, intronic 38, missense 60, nonsense 30, splice region 23, synonymous 6, utr 3. Merged into other: synonymous, utr. Model classes: canonical splice, intronic, missense, nonsense, other, splice region. `diagnostic_class_differs_from_cases_tsv` and `case_variants_without_mane_cadd_row` are 0 in every run.

## Invariance

`tests/test_invariance.py` reruns an analysis with one score dropped and compares every row that does not involve it (`logs/cb_invariance.log`): replication main without GPN-Star-M, 2,703 rows compared, 0 problems; without AVI, 2,357 rows, 0 problems; replication evee without GPN-MSA, 725 rows, 0 problems.
