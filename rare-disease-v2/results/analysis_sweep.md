# Frequency sweep (exploratory)

PROTOCOL.md section 13. Added on 6 November 2026 after every other result had been seen, so it is exploratory, not confirmatory. Made by `case_benchmark.py --study <study> --analysis sweep` from the 5%-cap build in `data/af5/`: `results/<study>_sweep.json`.

At each threshold: background grpmax (family left out) below it, diagnostics full gnomAD grpmax below it (a case whose diagnostic fails is dropped), 1000 Genomes AF below the larger of the threshold and 1%, then the inheritance filter. Intervals are 1,000 resamples of people; the 0.1% and 1% rows use the main and bg1pct draws, so they reproduce those results.

## Replication study, 500 cases built

Top 1, % (95% interval). Combined: lookups combined (locked, fitted on the 200 under the same threshold). Gains in percentage points with paired 95% intervals. Evo 2: not scored above 1%.

| Threshold | Cases | Median list | AVI | CADD | SpliceAI | GPN-Star-V | Evo 2 | Rule-based | Combined | Rule + AVI minus Rule | Combined + AVI minus combined |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01% | 464 (36 lost) | 4 | 92.7 (87.8 to 96.7) | 84.0 (76.3 to 90.2) | 72.6 (66.7 to 78.6) | 76.7 (68.0 to 84.4) | 78.7 (71.0 to 86.1) | 95.9 (92.9 to 98.2) | 97.4 (95.5 to 98.9) | +1.7 (+0.2 to +4.0), p = 0.032 | 0.0 (-1.1 to +1.3), p = 1 |
| 0.05% | 493 (7 lost) | 11 | 84.6 (79.0 to 89.4) | 65.3 (57.8 to 72.4) | 60.2 (54.9 to 65.1) | 64.1 (56.7 to 71.0) | 64.9 (57.9 to 72.1) | 89.2 (84.8 to 93.1) | 92.1 (89.2 to 94.8) | +2.2 (0.0 to +4.9), p = 0.08 | +1.8 (+0.2 to +3.5), p = 0.04 |
| 0.1% | 500 | 19 | 83.6 (78.2 to 88.8) | 61.0 (53.4 to 68.8) | 57.8 (53.0 to 62.6) | 59.4 (51.8 to 66.0) | 58.8 (52.0 to 65.4) | 87.2 (82.8 to 91.0) | 90.6 (87.8 to 93.4) | +2.6 (+0.4 to +5.6), p = 0.034 | +1.8 (-0.4 to +4.0), p = 0.162 |
| 0.5% | 500 | 50.5 | 76.2 (70.6 to 81.8) | 48.8 (42.4 to 55.2) | 52.4 (47.2 to 57.0) | 50.0 (43.6 to 56.4) | 41.6 (35.4 to 47.8) | 77.8 (71.6 to 84.0) | 86.8 (83.0 to 91.0) | +6.0 (+2.2 to +10.4), p = 0.001 | -3.0 (-9.0 to +1.6), p = 0.236 |
| 1% | 500 | 69.5 | 71.4 (65.2 to 77.4) | 44.2 (38.2 to 50.8) | 49.6 (45.0 to 54.6) | 47.8 (42.0 to 54.0) | 36.6 (30.8 to 42.4) | 74.8 (68.0 to 81.0) | 82.0 (76.4 to 87.2) | +5.2 (+1.4 to +9.6), p = 0.002 | -1.2 (-6.6 to +3.2), p = 0.724 |
| 2% | 500 | 101 | 68.6 (62.0 to 75.0) | 39.0 (33.0 to 45.0) | 49.0 (44.0 to 54.2) | 45.4 (39.0 to 52.0) | not scored | 69.8 (61.8 to 76.4) | 82.8 (76.8 to 87.8) | +6.0 (+2.0 to +10.4), p = 0.001 | -4.4 (-10.0 to 0.0), p = 0.054 |
| 5% | 500 | 177 | 63.2 (56.2 to 69.6) | 36.2 (30.2 to 42.2) | 47.2 (42.4 to 52.4) | 39.4 (33.6 to 45.6) | not scored | 61.4 (54.6 to 68.2) | 73.6 (67.6 to 79.4) | +6.4 (+2.8 to +10.4), p = 0.001 | -1.2 (-7.0 to +4.2), p = 0.784 |

Other single scores, top 1 %:

| Threshold | REVEL | AlphaMissense | GPN-MSA | GPN-Star-M |
|---|---|---|---|---|
| 0.01% | 48.7 (40.0 to 57.5) | 44.6 (35.6 to 54.1) | 74.1 (65.2 to 82.2) | 75.0 (66.5 to 83.4) |
| 0.05% | 29.8 (25.8 to 34.0) | 24.5 (21.0 to 28.3) | 56.4 (47.7 to 63.7) | 59.6 (51.1 to 67.4) |
| 0.1% | 28.6 (24.4 to 32.4) | 24.2 (20.6 to 28.0) | 51.2 (43.8 to 58.2) | 55.8 (48.0 to 63.4) |
| 0.5% | 27.0 (23.0 to 31.0) | 23.6 (19.8 to 27.2) | 38.2 (31.6 to 44.8) | 45.2 (39.0 to 51.4) |
| 1% | 25.6 (21.8 to 29.4) | 22.2 (19.0 to 25.6) | 36.0 (29.4 to 42.4) | 40.2 (33.8 to 46.6) |
| 2% | 25.6 (21.4 to 29.4) | 21.8 (18.2 to 25.0) | 32.0 (26.0 to 38.4) | 38.0 (32.0 to 45.2) |
| 5% | 23.6 (19.8 to 27.6) | 20.8 (17.4 to 24.2) | 27.0 (22.2 to 32.0) | 33.8 (28.2 to 39.4) |

AVI alone by whether the diagnostic is in full gnomAD (AC above 0), top 1 %, and the gains in each (paired; absent minus found on shared draws):

| Threshold | Found: cases | Found: AVI | Absent: cases | Absent: AVI | Rule + AVI gain, found | Rule + AVI gain, absent | Absent minus found |
|---|---|---|---|---|---|---|---|
| 0.01% | 78 | 88.5 (79.5 to 96.1) | 386 | 93.5 (88.6 to 97.4) | +2.6 (0.0 to +6.3) | +1.6 (0.0 to +3.4) | -1.0 (-4.6 to +1.9) |
| 0.05% | 107 | 78.5 (68.9 to 87.0) | 386 | 86.3 (80.8 to 92.0) | +1.9 (-1.9 to +7.0) | +2.3 (0.0 to +5.2) | +0.5 (-4.5 to +4.6) |
| 0.1% | 114 | 75.4 (65.6 to 84.2) | 386 | 86.0 (80.0 to 91.2) | +3.5 (0.0 to +8.2) | +2.3 (0.0 to +5.2) | -1.2 (-5.3 to +2.8) |
| 0.5% | 114 | 66.7 (57.0 to 75.4) | 386 | 79.0 (72.6 to 84.5) | +8.8 (+2.6 to +16.1) | +5.2 (+1.3 to +9.7) | -3.6 (-9.5 to +2.1) |
| 1% | 114 | 61.4 (51.7 to 70.2) | 386 | 74.4 (67.6 to 80.7) | +8.8 (+2.6 to +16.2) | +4.2 (0.0 to +8.4) | -4.6 (-10.1 to +0.8) |
| 2% | 114 | 58.8 (48.7 to 68.1) | 386 | 71.5 (64.5 to 77.8) | +9.7 (+3.5 to +16.3) | +4.9 (+1.0 to +8.8) | -4.7 (-10.2 to +0.5) |
| 5% | 114 | 57.0 (45.3 to 67.8) | 386 | 65.0 (57.8 to 72.2) | +8.8 (+2.8 to +16.0) | +5.7 (+2.3 to +9.8) | -3.1 (-8.3 to +1.9) |

Why: Spearman correlation of each score with background grpmax (distinct background SNVs of the ranked lists), and the mean number of background variants per case ranked at or above the spiked-in variant, by background grpmax band, next to the mean number in each band:

| Threshold | Background SNVs | AVI rho | CADD rho | In list, below 0.1% | In list, 0.1 to 1% | In list, 1 to 5% | Above by AVI, below 0.1% | Above by AVI, 0.1 to 1% | Above by AVI, 1 to 5% | Above by CADD, below 0.1% | Above by CADD, 0.1 to 1% | Above by CADD, 1 to 5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01% | 232 | -0.138 | -0.107 | 4.61 | 0.00 | 0.00 | 0.11 | 0.00 | 0.00 | 0.38 | 0.00 | 0.00 |
| 0.05% | 640 | -0.096 | -0.093 | 12.76 | 0.00 | 0.00 | 0.28 | 0.00 | 0.00 | 0.97 | 0.00 | 0.00 |
| 0.1% | 1002 | -0.095 | -0.101 | 19.95 | 0.00 | 0.00 | 0.38 | 0.00 | 0.00 | 1.35 | 0.00 | 0.00 |
| 0.5% | 2711 | -0.070 | -0.088 | 21.19 | 32.80 | 0.00 | 0.40 | 0.46 | 0.00 | 1.43 | 1.58 | 0.00 |
| 1% | 4200 | -0.077 | -0.092 | 21.69 | 61.95 | 0.00 | 0.42 | 0.80 | 0.00 | 1.49 | 2.84 | 0.00 |
| 2% | 6849 | -0.076 | -0.100 | 22.61 | 65.78 | 48.21 | 0.46 | 0.86 | 0.40 | 1.57 | 3.00 | 1.81 |
| 5% | 13581 | -0.071 | -0.095 | 24.73 | 71.30 | 175.08 | 0.55 | 0.90 | 1.67 | 1.76 | 3.24 | 6.11 |

Reproduction of the v6 results by the 0.1% and 1% rows (every value both compute: cases, list length, coverage, ties, every score and contrast in all cases and in the gnomAD strata, and the gnomAD split):

- 0.1% against `replication_main.json`: identical, 487 values compared, 0 differ, 0 missing from the reference.
- 1% against `replication_bg1pct.json`: identical, 487 values compared, 0 differ, 0 missing from the reference.

## SNV study, 200 cases built

Top 1, % (95% interval). Combined: lookups combined (cross-validated). Gains in percentage points with paired 95% intervals. Evo 2: not scored above 1%.

| Threshold | Cases | Median list | AVI | CADD | SpliceAI | GPN-Star-V | Evo 2 | Rule-based | Combined | Rule + AVI minus Rule | Combined + AVI minus combined |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01% | 188 (12 lost) | 8 | 86.7 (79.7 to 92.5) | 74.5 (65.8 to 83.3) | 63.3 (55.0 to 71.5) | 73.9 (65.0 to 82.4) | 72.3 (63.4 to 80.7) | 86.7 (78.1 to 93.4) | 91.0 (84.5 to 96.3) | +3.7 (+0.5 to +7.5), p = 0.048 | 0.0 (-2.2 to +2.7), p = 1 |
| 0.05% | 198 (2 lost) | 17 | 82.8 (75.9 to 89.3) | 57.1 (48.7 to 65.7) | 56.1 (48.0 to 62.9) | 58.6 (50.8 to 66.2) | 61.1 (53.1 to 69.2) | 84.9 (77.4 to 91.5) | 83.3 (74.4 to 91.0) | +4.5 (+1.0 to +8.1), p = 0.01 | +3.5 (-1.6 to +8.2), p = 0.236 |
| 0.1% | 200 | 26 | 79.0 (71.0 to 86.0) | 51.0 (43.0 to 59.5) | 53.5 (46.0 to 59.5) | 52.5 (44.0 to 60.5) | 51.5 (43.5 to 59.5) | 81.0 (73.5 to 88.0) | 84.0 (75.5 to 91.5) | +4.5 (+0.5 to +8.0), p = 0.028 | +4.0 (-2.0 to +9.5), p = 0.178 |
| 0.5% | 200 | 58.5 | 72.0 (63.5 to 80.0) | 41.5 (34.0 to 49.5) | 48.5 (42.0 to 55.0) | 46.0 (39.0 to 53.0) | 44.0 (36.5 to 51.5) | 74.5 (66.0 to 82.5) | 75.5 (66.5 to 84.5) | +3.5 (-0.5 to +7.5), p = 0.082 | +7.0 (-1.5 to +15.0), p = 0.12 |
| 1% | 200 | 78.5 | 70.0 (61.5 to 78.0) | 41.0 (33.5 to 48.5) | 47.0 (40.0 to 54.5) | 45.0 (38.0 to 52.0) | 39.0 (32.5 to 46.0) | 71.5 (62.5 to 79.5) | 77.0 (68.0 to 85.0) | +4.0 (+0.5 to +8.0), p = 0.038 | +1.0 (-7.5 to +9.5), p = 0.85 |
| 2% | 200 | 114 | 65.0 (56.5 to 72.5) | 38.0 (29.0 to 46.5) | 45.0 (37.5 to 52.5) | 42.0 (34.5 to 49.5) | not scored | 67.0 (58.5 to 75.5) | 72.0 (63.0 to 80.5) | +3.5 (-1.5 to +8.0), p = 0.2 | +6.0 (-3.0 to +14.5), p = 0.21 |
| 5% | 200 | 199.5 | 62.0 (54.0 to 70.0) | 38.0 (30.5 to 46.5) | 44.0 (36.5 to 51.5) | 36.0 (29.0 to 43.5) | not scored | 65.0 (56.5 to 74.0) | 66.5 (57.5 to 76.5) | +2.5 (-2.0 to +7.0), p = 0.298 | +4.0 (-6.0 to +13.5), p = 0.46 |

Other single scores, top 1 %:

| Threshold | REVEL | AlphaMissense | GPN-MSA | GPN-Star-M |
|---|---|---|---|---|
| 0.01% | 38.8 (30.1 to 48.1) | 36.7 (27.7 to 46.1) | 69.7 (59.7 to 79.2) | 70.7 (60.9 to 79.9) |
| 0.05% | 28.8 (23.5 to 33.8) | 25.8 (21.5 to 30.1) | 52.0 (42.1 to 61.3) | 57.6 (49.2 to 66.0) |
| 0.1% | 27.5 (23.0 to 32.5) | 23.5 (19.0 to 28.0) | 45.5 (36.0 to 55.0) | 50.5 (42.5 to 58.0) |
| 0.5% | 26.0 (21.0 to 31.0) | 22.5 (18.0 to 27.0) | 37.0 (28.0 to 46.0) | 41.0 (34.0 to 48.5) |
| 1% | 25.5 (20.5 to 30.5) | 22.5 (18.0 to 27.5) | 34.5 (27.5 to 42.5) | 38.0 (30.0 to 46.0) |
| 2% | 24.0 (19.0 to 29.0) | 21.5 (16.5 to 26.0) | 28.5 (21.5 to 36.5) | 34.5 (26.0 to 42.5) |
| 5% | 23.0 (18.0 to 28.5) | 21.5 (17.0 to 26.0) | 26.5 (19.0 to 34.0) | 28.5 (20.5 to 36.5) |

AVI alone by whether the diagnostic is in full gnomAD (AC above 0), top 1 %, and the gains in each (paired; absent minus found on shared draws):

| Threshold | Found: cases | Found: AVI | Absent: cases | Absent: AVI | Rule + AVI gain, found | Rule + AVI gain, absent | Absent minus found |
|---|---|---|---|---|---|---|---|
| 0.01% | 40 | 85.0 (74.4 to 94.7) | 148 | 87.2 (78.1 to 93.8) | 0.0 (0.0 to 0.0) | +4.7 (+0.7 to +9.7) | +4.7 (+0.7 to +9.7) |
| 0.05% | 50 | 76.0 (63.0 to 87.8) | 148 | 85.1 (77.3 to 92.4) | 0.0 (0.0 to 0.0) | +6.1 (+1.4 to +10.9) | +6.1 (+1.4 to +10.9) |
| 0.1% | 52 | 73.1 (60.3 to 86.0) | 148 | 81.1 (71.6 to 89.8) | 0.0 (-6.7 to +5.0) | +6.1 (+1.9 to +11.5) | +6.1 (-0.6 to +13.8) |
| 0.5% | 52 | 63.5 (50.0 to 77.1) | 148 | 75.0 (65.6 to 84.2) | -1.9 (-9.1 to +3.9) | +5.4 (+0.7 to +10.3) | +7.3 (-0.4 to +15.1) |
| 1% | 52 | 61.5 (47.9 to 74.0) | 148 | 73.0 (63.5 to 81.2) | 0.0 (-4.7 to +5.3) | +5.4 (+0.7 to +10.4) | +5.4 (-1.2 to +12.7) |
| 2% | 52 | 57.7 (44.4 to 71.7) | 148 | 67.6 (58.7 to 76.3) | -3.9 (-11.1 to +3.3) | +6.1 (+1.2 to +11.7) | +9.9 (+1.6 to +18.1) |
| 5% | 52 | 57.7 (42.9 to 72.9) | 148 | 63.5 (54.5 to 72.7) | -3.9 (-11.3 to +3.6) | +4.7 (0.0 to +9.1) | +8.6 (+0.8 to +16.9) |

Why: Spearman correlation of each score with background grpmax (distinct background SNVs of the ranked lists), and the mean number of background variants per case ranked at or above the spiked-in variant, by background grpmax band, next to the mean number in each band:

| Threshold | Background SNVs | AVI rho | CADD rho | In list, below 0.1% | In list, 0.1 to 1% | In list, 1 to 5% | Above by AVI, below 0.1% | Above by AVI, 0.1 to 1% | Above by AVI, 1 to 5% | Above by CADD, below 0.1% | Above by CADD, 0.1 to 1% | Above by CADD, 1 to 5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01% | 405 | -0.075 | -0.030 | 8.05 | 0.00 | 0.00 | 0.29 | 0.00 | 0.00 | 0.83 | 0.00 | 0.00 |
| 0.05% | 884 | -0.103 | -0.054 | 17.67 | 0.00 | 0.00 | 0.48 | 0.00 | 0.00 | 1.55 | 0.00 | 0.00 |
| 0.1% | 1301 | -0.041 | -0.019 | 26.00 | 0.00 | 0.00 | 0.73 | 0.00 | 0.00 | 2.31 | 0.00 | 0.00 |
| 0.5% | 3123 | -0.087 | -0.099 | 27.41 | 34.98 | 0.00 | 0.80 | 0.99 | 0.00 | 2.40 | 2.17 | 0.00 |
| 1% | 4844 | -0.078 | -0.110 | 28.39 | 68.42 | 0.00 | 0.85 | 2.12 | 0.00 | 2.48 | 3.89 | 0.00 |
| 2% | 7730 | -0.079 | -0.089 | 29.86 | 72.06 | 52.54 | 0.90 | 2.23 | 1.62 | 2.58 | 4.08 | 2.87 |
| 5% | 14980 | -0.094 | -0.105 | 32.02 | 77.53 | 189.86 | 0.98 | 2.40 | 6.15 | 2.75 | 4.42 | 8.71 |

Reproduction of the v6 results by the 0.1% and 1% rows (every value both compute: cases, list length, coverage, ties, every score and contrast in all cases and in the gnomAD strata, and the gnomAD split):

- 0.1% against `snv_main.json`: identical, 414 values compared, 0 differ, 0 missing from the reference.
- 1% against `snv_bg1pct.json`: identical, 414 values compared, 0 differ, 0 missing from the reference.
