# Replication: does GPN-MSA's top-1 gain hold on new cases?

Raised in review of the talk (issue 19 in the private decks repo). The original 200-case SNV benchmark (`results.json`) put the diagnostic variant at top 1 in 59.5% of
cases with the four lookups combined (CADD v1.7, SpliceAI, REVEL, AlphaMissense, plus consequence class) and 65.5%
with GPN-MSA added: +6.0 points, 95% CI +2.0 to +10.5. Top 5 and top 10 did not change. This is a test of that
result on 500 new cases.

## Method

- **Diagnostic variants.** 500 ClinVar 2-star P/LP SNVs in the same 1,512-gene intellectual-disability panel,
  none at a position used as a diagnostic in the original 200. Same consequence mix scaled to 500: 150 missense,
  100 canonical splice, 75 nonsense, 175 synonymous, UTR or intronic. Seed 20261023
  (`build_cases.py replication`).
- **Backgrounds.** 50 new 1000 Genomes high-coverage samples, drawn with seed 20261023 from the 2,504 unrelated
  samples, leaving out the original 50 and anyone in their families (2,425 eligible). The original 50 included 16
  of the 698 related samples; the new 50 include none. Superpopulations: EAS 13, EUR 12, AMR 9, AFR 8, SAS 8
  (original: AFR 12, EUR 11, EAS 11, SAS 10, AMR 6). Genotypes fetched by `fetch_kgp.sh` with the new sample file,
  except chrX and chr15. chrX was too slow as one remote stream and chr15 came back truncated at 25 Mb while still
  printing "done" (the script hides bcftools errors and the pipeline's exit status is the last command's). Both
  were refetched as 8 parallel region chunks and merged with `bcftools concat | bcftools norm -d exact`. Every
  chromosome was checked: 50 samples, last record within the last panel region, and 85 to 95% of the original
  record count. The original 50-sample files pass the same end-position check.
  Each sample carries 10 cases. Fewer African backgrounds give fewer rare variants: median 208.5 per case
  (original 234.5).
- **Scores.** CADD v1.7 (with SpliceAI), REVEL, AlphaMissense, GPN-MSA, GPN-Star vertebrate (V) and mammal (M), by
  `lookup_scores.py --study replication --reuse`, which copies rows for the 2,699 of 11,708 variants already
  scored in the original study. Evo 2 was not rerun.
- **Coverage of diagnostic variants.** CADD, GPN-MSA, GPN-Star V and M: 500/500 (100%). SpliceAI 475, REVEL 159,
  AlphaMissense 147 (of 150 missense).
- **Locked models.** "Lookups combined", "+ GPN-MSA", "+ GPN-Star-V" and "+ GPN-Star-M" are fitted once on all 200
  original cases, with the same features and logistic regression as `case_benchmark.py`, then applied unchanged to
  the 500 new cases.
- **Cross-validated.** The original analysis repeated on the 500 new cases alone: 5-fold cross-validation grouped
  by background sample.
- **Intervals.** 95% bootstrap intervals, 1,000 resamples of the 50 background samples (each with its 10 cases),
  as in the original. With 50 clusters the intervals are wider than 500 independent cases would give.

`case_benchmark.py --study replication` writes `results_replication.json`. Rerunning the original analysis from
this copy reproduces every lookup and GPN row of `results.json` exactly (scikit-learn 1.9.1); the Evo 2 rows
differ by one case.

## Results (500 cases, % with the diagnostic variant at or above rank k)

Locked models, fitted on the original 200 cases:

| Model | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| Lookups combined | 61.0 [54.2, 67.8] | 91.0 [88.8, 93.2] | 94.0 [92.0, 96.0] |
| Lookups + GPN-MSA | 62.2 [54.4, 69.4] | 91.6 [89.4, 93.8] | 94.2 [92.0, 96.0] |
| Lookups + GPN-Star-V | 64.0 [57.0, 70.0] | 91.6 [89.6, 93.4] | 94.4 [92.6, 96.0] |
| Lookups + GPN-Star-M | 63.6 [56.4, 70.2] | 91.8 [89.4, 94.0] | 96.0 [94.4, 97.4] |

Cross-validated on the 500 new cases:

| Model | Top 1 | Top 5 | Top 10 |
|---|---|---|---|
| Lookups combined | 58.2 [50.6, 65.2] | 92.6 [90.4, 94.6] | 95.8 [94.0, 97.4] |
| Lookups + GPN-MSA | 60.6 [53.4, 67.0] | 93.6 [91.6, 95.4] | 96.0 [94.4, 97.4] |
| Lookups + GPN-Star-V | 63.8 [57.0, 70.2] | 93.8 [92.4, 95.4] | 96.8 [95.2, 98.2] |
| Lookups + GPN-Star-M | 65.0 [59.2, 70.8] | 94.0 [92.0, 95.8] | 96.4 [94.6, 98.0] |

Gain over lookups combined, same cases (points):

| Added score | Original, top 1 | Locked, top 1 | CV, top 1 | Locked, top 5 | Locked, top 10 |
|---|---|---|---|---|---|
| GPN-MSA | +6.0 [+2.0, +10.5] | +1.2 [-2.8, +5.2] | +2.4 [-0.6, +5.8] | +0.6 [-1.2, +2.4] | +0.2 [-1.0, +1.6] |
| GPN-Star-V | +4.0 [-2.0, +10.0] | +3.0 [-1.2, +7.8] | +5.6 [+1.8, +9.6] | +0.6 [-1.0, +2.2] | +0.4 [-1.0, +1.8] |
| GPN-Star-M | +1.5 [-2.0, +5.5] | +2.6 [-2.2, +6.8] | +6.8 [+2.8, +11.8] | +0.8 [-0.6, +2.4] | +2.0 [+0.8, +3.2] |

One more observation. For missense diagnostics, the combined model is worse at top 1 than AlphaMissense alone
(replication CV: 35% against 67%; original: 43% against 58%). The combined model is not the best ranker
available at top 1, so some of what an added score "gains" there may be making up for that.

## Conclusion

The +6 point top-1 gain from adding GPN-MSA attenuates to the point of not replicating. On 500 new cases with new
backgrounds, +6.0 [+2.0, +10.5] becomes +1.2 [-2.8, +5.2] with the model locked from the original fit and +2.4
[-0.6, +5.8] when refitted on the new cases. Both intervals include zero, and the original estimate lies outside
the locked interval. The original estimate was the largest of four scores tested on 200 cases, so a
smaller value on new data is what one would expect. GPN-Star does a little better: +3 points at top 1 when locked
(interval includes zero) and +6 to +7 points when refitted on the larger set (interval excludes zero), so a small
top-1 gain from a conservation score is possible but not established; the refitted gains are the same "best of
several scores" pattern that produced the original GPN-MSA claim. At top 5 and top 10 none of the added scores
changes the result by more than 2 points. For the slides: the lookups combined with consequence do the work, and
the DNA foundation-model scores add little or nothing.
