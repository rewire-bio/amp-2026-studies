---
title: "Does Evo 2 help rank the loss-of-function variant in a simulated BRCA1 case?"
subtitle: "A case-level ranking test on saturation genome editing labels"
author: '`Tim Richardson \\ Rewire`{=latex}'
date: 2 November 2026
abstract: |
  We asked whether a zero-shot score from the Evo 2 7B genome language model helps put a loss-of-function variant
  near the top of a list of candidates. We used 3,893 BRCA1 single-nucleotide variants with saturation genome
  editing labels, of which 823 are loss of function. Each simulated case held one loss-of-function variant and 100
  functional or intermediate variants from other positions, giving 8,230 cases. We ranked each case by CADD v1.7,
  SpliceAI, REVEL v1.3, Evo 2 7B and cross-validated combinations of them, and report how often the
  loss-of-function variant was ranked first or in the top 5, with 95% intervals from resampling positions. On its
  own, Evo 2 put the variant in the top 5 in 48.1% of cases (44.1 to 52.6), against 54.9% (50.6 to 59.3) for CADD
  v1.7. Added to a combination of CADD and SpliceAI with separate weights for coding and noncoding variants, Evo 2
  raised the top-5 rate by 3.1 points (1.1 to 5.2). That combination's coding top-5 rate (51.3%) was below CADD
  alone (56.5%), though the two were not compared by a paired test. Against CADD alone on coding variants, CADD +
  Evo 2 gave +2.5 points (-1.3 to 6.8), no gain detected. For noncoding variants, no combination's point estimate
  reached SpliceAI alone (95.3%); this was not tested. Evo 2's apparent gain tended to be larger on weaker
  baselines, but the difference between gains was not tested. Seven of the 30 paired comparisons changed
  significance with 20 candidates per case instead of 100. All comparisons are exploratory.
---

# Introduction

A lab that sequences a patient finds many rare variants and has to decide which to review first. It ranks them
with established scores such as CADD [@schubach2024cadd], SpliceAI [@jaganathan2019predicting] and REVEL
[@ioannidis2016revel], usually inside a prioritisation tool such as Exomiser [@smedley2015next]. Genome language
models such as Evo 2 [@brixi2026genome] score a variant without training on labels, from how much the variant
changes the likelihood the model assigns to the surrounding DNA. Arc Institute showed this zero-shot score
separating loss-of-function from functional BRCA1 variants [@arcinstitute2026evo2], and a reproducibility study
of the same scores found no difference between Evo 2 7B and CADD v1.7 overall: AUROC 0.875 (0.859 to 0.891)
against 0.881 (0.862 to 0.897), difference -0.006 (-0.025 to 0.014) [@rewire2026evo2ops].

A variant-level AUROC does not say whether a score helps in a case. A lab looks at the top few candidates, so
what matters is whether the diagnostic variant is ranked first, or in the first few, among the others in the same
case. It also matters whether a new score adds anything to the scores the lab already uses.

This study asks that question on one gene where every variant has a functional label from a lab assay: the BRCA1
saturation genome editing data of Findlay et al. [@findlay2018accurate]. The labels do not depend on ClinVar or on
the predictors being tested. It was the first case-level test behind a talk at AMP 2026; a later rare-disease
benchmark with real backgrounds and ClinVar diagnoses replaced it [@rewire2026studies; @rewire2026rdv2].

The questions were:

1. How often does each score, alone, rank the loss-of-function variant first or in the top 5?
2. Does adding Evo 2 to CADD, or to CADD and SpliceAI combined, raise those rates?
3. Does the answer change by variant type (coding or noncoding), or with the way scores are combined?

# Methods

The script `brca1-case-ranking/case_ranking.py` is the full record. This section follows it.

## Variants and labels

The variants are the 3,893 BRCA1 single-nucleotide variants (SNVs) scored by saturation genome editing in 13
exons [@findlay2018accurate], as prepared in `data/derived/brca1_variants.tsv` of the Evo 2 cloud operations study
[@rewire2026evo2ops]. A variant is loss of function (LOF) if Findlay et al. class it `LOF`; the other 3,070
variants are classed functional (`FUNC`) or intermediate (`INT`) and are treated as non-diagnostic. Coordinates are
GRCh37. A variant is coding if its consequence is missense, nonsense or synonymous, and noncoding otherwise
(2,768 coding variants with 584 LOF; 1,125 noncoding with 239 LOF) [@rewire2026evo2ops]. 230 of the 239 noncoding
LOF variants are canonical splice or splice-region variants [@rewire2026evo2ops].

## Scores

Every score is oriented so that higher means more likely to be loss of function.

- **Evo 2 7B.** The `evo2_7b` checkpoint (Hugging Face `arcinstitute/evo2_7b`, revision `bda0089`, Apache-2.0
  [@arcinstitute2025evo2]) scored 8,192 bp windows of GRCh37 centred on each SNV, following Arc's BRCA1 notebook
  [@arcinstitute2026evo2]. The score is the reference window's mean log-likelihood minus the variant window's
  (the negated delta). The scores are those of the recorded run `b75c67aa` on a Google Cloud H100 80 GB, in the
  container `ghcr.io/rewire-bio/evo2-cloud-operations@sha256:22970720...` [@rewire2026evo2ops].
- **CADD v1.7.** The raw score from the GRCh37 `whole_genome_SNVs.tsv.gz` [@schubach2024cadd; @caddteam2026download].
- **SpliceAI.** Illumina's masked precomputed scores as annotated in CADD v1.7's `whole_genome_SNVs_inclAnno.tsv.gz`,
  BRCA1 rows only; the score is the largest of the four delta scores [@jaganathan2019predicting;
  @illumina2019spliceai].
- **REVEL v1.3.** The largest score over transcripts from `revel-v1.3_all_chromosomes.zip` [@ioannidis2016revel;
  @revelteam2026downloads]. REVEL is a missense score. It scores 2,144 variants: all 2,086 missense variants,
  50 of the 682 nonsense or synonymous variants and 8 of the 1,125 noncoding variants [@rewire2026evo2ops].

The CADD, SpliceAI and REVEL extracts for chr17:41197000-41277000 are fetched by
`followup/current_predictors.py --fetch` in the Evo 2 study [@rewire2026evo2ops]. That step writes the extracts'
SHA-256 sums to `evidence/current-predictors/sources.json`. An unchanged `git diff` of that file confirms the
extracts match the ones used here. CADD v1.7, SpliceAI and Evo 2 score all 3,893 variants; REVEL scores 2,144.

## Combinations

Three kinds of combination were built from the raw scores:

- **Logistic models.** A standardised logistic regression (scikit-learn `StandardScaler` then
  `LogisticRegression(max_iter=1000)`, default L2 penalty) on CADD + SpliceAI, CADD + Evo 2, and CADD + SpliceAI +
  Evo 2 [@scikitlearn2026linear].
- **Type-aware logistic models.** The same, with each score also multiplied by a coding indicator, and the
  indicator itself, so each score carries a different weight for coding and noncoding variants. Built for
  CADD + SpliceAI and CADD + SpliceAI + Evo 2.
- **Average of percentiles.** Each score's percentile rank over all 3,893 variants, averaged, for CADD + SpliceAI and
  CADD + SpliceAI + Evo 2. This uses no labels.

Every logistic model's score is out of fold: 5-fold cross-validation with `GroupKFold` grouped by genomic position
[@scikitlearn2026groupkfold], so no variant is scored by a model that saw its label or any other allele at its
position.

## Simulated cases

Each of the 823 LOF variants is the diagnostic variant in 10 cases. Each case adds K candidates drawn without
replacement from the 3,070 non-LOF variants at other positions than the diagnostic variant. K = 100 is the main
analysis (8,230 cases of 101 variants); K = 20 is a second case size (8,230 cases of 21 variants). Draws use a
NumPy generator seeded with 20261016 + K.

For each score, a competitor is a candidate that scores equal to or higher than the diagnostic variant. Ties count
against the score. A missing score ranks lowest: a candidate without a score never outranks, and a diagnostic
variant without a score is outranked by every candidate with one. The diagnostic variant is ranked first if it has
no competitors and is in the top 5 if it has at most 4. We report both rates for all cases, and for cases whose
diagnostic variant is coding (5,840) or noncoding (2,390).

## Statistics

Intervals are 95% percentile intervals from 1,000 bootstrap resamples [@efron1979bootstrap] of diagnostic variants
grouped by position: each resample draws positions with replacement and keeps every case whose diagnostic variant is
at a drawn position. Within a stratum, only that stratum's positions are drawn. Intervals for single rates use
separate resamples for each score. Paired differences ("adding Evo 2") compare two scores on the same cases and use
the same resamples for both. A difference is called significant at the 5% level when its 95% interval excludes
zero. There is no correction for multiple comparisons. Bootstrap resamples come from the same generator as the case
draws, after them. For each case size, the code works through the strata in turn (all, coding, noncoding) and, within
each, draws every single-score interval and then every paired interval, in the order the scores and pairs appear in
`case_ranking.py`. Adding or reordering a score therefore changes later intervals.

Five paired comparisons are coded: CADD + Evo 2 against CADD; CADD + SpliceAI + Evo 2 against CADD + SpliceAI;
CADD + SpliceAI + Evo 2 against CADD (which adds SpliceAI and Evo 2 together); and, for the type-aware and average-of-percentile combinations, the
version with Evo 2 against the one without. Each is computed for both rates, three strata and two case sizes. No
analysis plan was written before the results were seen, and the combination methods were developed one after
another as earlier results came in (the average of percentiles was tried first and the type-aware model last).
All comparisons are therefore exploratory. The type-aware comparison on all cases at K = 100 is the headline
because it has the strongest baseline, not because it was fixed in advance. Single scores were not compared with
each other by paired tests.

## Reproducing the run

```
uv run --with pandas --with numpy --with scikit-learn \
  python brca1-case-ranking/case_ranking.py --study ../evo2-cloud-operations
```

Clone `rewire-bio/evo2-cloud-operations` next to this repository and check out commit `f8ae11c`. In that checkout
run `make inputs` (writes `data/derived/brca1_variants.tsv`) and
`uv run --frozen python followup/current_predictors.py --fetch` (needs tabix, curl, unzip and network access).
Confirm that `git diff evidence/current-predictors/sources.json` is empty. Then run the command above from this
repository's root. The command does not pin package versions. Rerun on 2 November 2026 on an Apple M4 laptop (CPU only) with Python
3.11.13, pandas 3.1.0, NumPy 2.4.6 and scikit-learn 1.9.1, it took about 7 seconds and wrote a `results.json`
byte-identical to the committed one. The Evo 2 scores needed a GPU: the H100 run scored 2.93 sequences per
second at a peak of 16.8 GiB of GPU memory [@rewire2026evo2ops]. The charts are drawn from `results.json` by
`paper/charts.py` (`bash brca1-case-ranking/paper/make_figures.sh`).

The tables in the study's README come from a run made before the average-of-percentiles combination was added.
Adding it changed the bootstrap draws (the simulated cases are unchanged), so some README interval bounds differ
from `results.json` by up to 0.6 points; no point estimate changed. This paper uses `results.json` throughout.

![Study design. Each LOF variant is the diagnostic variant in 10 simulated cases; the other candidates are
functional or intermediate variants at other positions.](figures/diagram_design.pdf){width=50%}

![How one case is scored. The rule for ties and missing scores is applied the same way to every score and
combination.](figures/diagram_test_case.pdf){width=62%}

![Code and data. The inputs come from the Evo 2 cloud operations study; everything after `case_ranking.py` is
in this repository.](figures/diagram_pipeline.pdf){width=90%}

# Results

## Single scores and combinations

With 100 candidates per case, CADD v1.7 was the best single score: it ranked the LOF variant first in 30.2% of
cases (26.9 to 33.7) and in the top 5 in 54.9% (50.6 to 59.3) (Table 1, Figure 4). Evo 2 ranked it first in 18.1%
(15.0 to 20.9) and in the top 5 in 48.1% (44.1 to 52.6). The ranked-first intervals of CADD and Evo 2 do not
overlap; the top-5 intervals do, and no paired test between them was run. SpliceAI and REVEL were lower overall
because each is weak outside its own variant type: SpliceAI on coding variants, and REVEL on noncoding variants.
REVEL scores 8 of the 1,125 noncoding variants, so it almost never ranks a noncoding LOF variant highly: it ranked
one first in no case, and in the top 5 in 0.0% of noncoding cases at K = 100 and 0.2% (0.04 to 0.35) at K = 20.

The best combination was the type-aware CADD + SpliceAI + Evo 2 model: 31.3% ranked first (28.1 to 34.5) and
65.0% in the top 5 (61.2 to 69.0). The average of percentiles of CADD and SpliceAI was below CADD alone in the top 5
(41.0% against 54.9%).

![Share of cases with the LOF variant ranked first (left) and in the top 5 (right), 100 candidates per case, with
95% intervals. Blue: no Evo 2. Orange: includes Evo 2. C is CADD v1.7 and S is SpliceAI.](figures/scores_k100.pdf){width=100%}

: Share of the 8,230 cases (K = 100) in which the LOF variant is ranked first or in the top 5, with 95% intervals,
and the median number of competitors. Source: `results.json`, `cases.k100.all`.

| Score | Ranked first, % (95% CI) | In top 5, % (95% CI) | Median competitors |
|----------------------------------------|----------------------:|----------------------:|--------------:|
| REVEL v1.3 | 13.7 (10.9 to 16.4) | 30.4 (26.6 to 34.7) | 24 |
| SpliceAI | 17.2 (14.2 to 20.4) | 37.2 (32.4 to 41.8) | 21 |
| Evo 2 7B | 18.1 (15.0 to 20.9) | 48.1 (44.1 to 52.6) | 5 |
| CADD v1.7 | 30.2 (26.9 to 33.7) | 54.9 (50.6 to 59.3) | 3 |
| Average percentile CADD + SpliceAI | 21.5 (17.9 to 25.3) | 41.0 (36.6 to 45.5) | 7 |
| Average percentile CADD + SpliceAI + Evo 2 | 21.9 (18.4 to 25.4) | 47.0 (42.6 to 51.4) | 5 |
| Model CADD + SpliceAI | 27.7 (24.4 to 31.2) | 57.4 (53.7 to 61.3) | 3 |
| Model CADD + Evo 2 | 27.6 (24.5 to 31.0) | 59.1 (55.2 to 63.3) | 3 |
| Model CADD + SpliceAI + Evo 2 | 29.8 (26.7 to 33.1) | 62.1 (58.3 to 65.9) | 2 |
| Type-aware CADD + SpliceAI | 30.6 (27.5 to 33.7) | 62.0 (58.2 to 65.7) | 2 |
| Type-aware CADD + SpliceAI + Evo 2 | 31.3 (28.1 to 34.5) | 65.0 (61.2 to 69.0) | 2 |

## Adding Evo 2

Table 2 and Figure 5 give the paired change from adding Evo 2 to each baseline. On all cases, adding Evo 2 raised
the top-5 rate against every baseline, and each of those gains is significant at the 5% level. On all cases the gain
was +6.0 points (3.9 to 7.9) on the average of percentiles and +3.1 (1.1 to 5.2) on the type-aware model. The
difference between gains was not tested, and on coding cases the order reversed (+4.3 against +4.8). For ranked first, adding Evo 2 to the type-aware model gave +0.8
points (-0.6 to 2.1), not significant. Adding Evo 2 to CADD gave -2.7 points (-5.7 to 0.4) for ranked first, not
significant, and +4.2 points (0.9 to 7.5) for the top 5.

![Paired change from adding Evo 2 to each baseline, by type of diagnostic variant, 100 candidates per case, with
95% intervals. Orange: interval above zero. Blue: interval below zero. Grey: interval includes zero.](figures/gains_k100.pdf){width=100%}

: Paired change in percentage points from adding Evo 2, K = 100, with 95% intervals. The rows "CADD (adds S + E)"
compare CADD + SpliceAI + Evo 2 with CADD alone, so they add SpliceAI and Evo 2 together and are not in Figure 5. Significant: the interval excludes zero (5% level,
no multiplicity correction; all exploratory). Source: `results.json`, `cases.k100.*.paired`.

| Cases | Baseline | Ranked first | Sig. | In top 5 | Sig. |
|----------|-------------------------|-----------------------:|----------|-----------------------:|----------|
| All | CADD alone | -2.7 (-5.7 to 0.4) | No | +4.2 (0.9 to 7.5) | Yes |
| All | CADD (adds S + E) | -0.4 (-2.9 to 2.2) | No | +7.2 (3.7 to 11.1) | Yes |
| All | Avg percentile C + S | +0.4 (-1.6 to 2.2) | No | +6.0 (3.9 to 7.9) | Yes |
| All | Model C + S | +2.1 (0.5 to 3.6) | Yes | +4.8 (2.8 to 7.1) | Yes |
| All | Type-aware C + S | +0.8 (-0.6 to 2.1) | No | +3.1 (1.1 to 5.2) | Yes |
| Coding | CADD alone | -2.2 (-5.2 to 0.8) | No | +2.5 (-1.3 to 6.8) | No |
| Coding | CADD (adds S + E) | -8.4 (-11.0 to -6.1) | Yes, loss | -2.0 (-5.8 to 1.7) | No |
| Coding | Avg percentile C + S | +2.9 (1.4 to 4.2) | Yes | +4.3 (2.2 to 6.3) | Yes |
| Coding | Model C + S | +3.1 (1.7 to 4.5) | Yes | +6.3 (3.3 to 9.0) | Yes |
| Coding | Type-aware C + S | +1.6 (0.1 to 3.1) | Yes | +4.8 (2.0 to 7.7) | Yes |
| Noncoding | CADD alone | -3.8 (-10.8 to 4.0) | No | +8.1 (2.7 to 14.0) | Yes |
| Noncoding | CADD (adds S + E) | +18.9 (13.8 to 24.0) | Yes | +29.8 (22.1 to 38.0) | Yes |
| Noncoding | Avg percentile C + S | -5.8 (-11.2 to -0.8) | Yes, loss | +10.1 (5.6 to 14.9) | Yes |
| Noncoding | Model C + S | -0.5 (-4.1 to 3.5) | No | +0.9 (-2.9 to 5.7) | No |
| Noncoding | Type-aware C + S | -1.3 (-3.7 to 1.4) | No | -1.1 (-2.0 to -0.2) | Yes, loss |

## By variant type

Splitting by the diagnostic variant's type shows where the gains come from (Table 3, Figure 6). For coding cases,
CADD v1.7 alone put the LOF variant in the top 5 in 56.5% (51.7 to 60.9). The type-aware CADD + SpliceAI model's
point estimate was lower, 51.3% (46.9 to 56.1); the two were not compared by a paired test. Adding Evo 2 to that model
brought it to 56.1% (52.0 to 60.6), close to CADD alone. Against CADD alone, the CADD + Evo 2 model gave +2.5 points
(-1.3 to 6.8): no gain detected. For noncoding cases, SpliceAI alone reached 95.3% (91.5 to 98.1); no combination's
point estimate reached it, though this was not tested. Adding Evo 2 to the type-aware model cost 1.1 points (-2.0 to
-0.2) at K = 100, but this loss was not seen at K = 20, where the same comparison gave +0.6 (0.04 to 1.34).

![Share of cases with the LOF variant in the top 5, by type of diagnostic variant, 100 candidates per case, with
95% intervals.](figures/types_k100.pdf){width=100%}

: Share of cases with the LOF variant in the top 5 by type of diagnostic variant, K = 100, with 95% intervals.
Source: `results.json`, `cases.k100.coding` and `cases.k100.noncoding`.

| Score | Coding (5,840 cases) | Noncoding (2,390 cases) |
|--------------------------------------------|----------------------------:|----------------------------:|
| CADD v1.7 | 56.5 (51.7 to 60.9) | 51.0 (40.7 to 60.1) |
| SpliceAI | 13.4 (10.5 to 16.9) | 95.3 (91.5 to 98.1) |
| REVEL v1.3 | 42.8 (38.2 to 47.4) | 0.0 (0.0 to 0.0) |
| Evo 2 7B | 45.0 (40.1 to 49.3) | 55.9 (46.7 to 64.6) |
| Model CADD + Evo 2 | 59.0 (54.7 to 63.7) | 59.1 (50.8 to 67.8) |
| Type-aware CADD + SpliceAI | 51.3 (46.9 to 56.1) | 88.0 (82.4 to 92.9) |
| Type-aware CADD + SpliceAI + Evo 2 | 56.1 (52.0 to 60.6) | 87.0 (80.7 to 91.7) |

SpliceAI is a splicing model and separates coding LOF variants poorly (variant-level AUROC 0.588, 0.559 to 0.618
[@rewire2026evo2ops]). With ties counted against the score, the median coding case had 45 SpliceAI competitors.

## Twenty candidates per case

With 20 candidates the rates are higher (Table 4). Seven of the 30 paired comparisons change significance between
the two case sizes, and two of them change sign. On noncoding cases, adding Evo 2 to the type-aware model lowered the
top-5 rate at K = 100 (-1.1 points, -2.0 to -0.2) but raised it at K = 20 (+0.6, 0.04 to 1.34), so neither result is
stable. Adding Evo 2 to the average of percentiles on noncoding cases changed ranked first by -5.8 points (-11.2 to
-0.8) at K = 100 and by +4.2 (0.4 to 8.1) at K = 20; its top-5 gain of +10.1 at K = 100 shrank to +1.4 (-0.5 to 3.3).

On all cases, three ranked-first comparisons became significant at K = 20: adding Evo 2 to the type-aware model, +2.6
points (1.1 to 4.2); adding it to the average of percentiles, +5.1 (3.5 to 6.6); and CADD + SpliceAI + Evo 2 against
CADD alone, +5.2 (2.0 to 8.2). Adding Evo 2 to the type-aware model also raised the top 5 by 0.9 points (0.2 to 1.7),
as at K = 100. For coding cases, the CADD + Evo 2 model raised the top-5 rate over CADD alone by 1.5 points (0.1 to
2.9), where K = 100 detected no gain. The gain from adding Evo 2 on all cases, top 5, was +11.9 points on the average
of percentiles (baseline 74.6%), +5.8 on CADD alone (79.6%), +1.4 on the logistic CADD + SpliceAI model (88.7%) and
+0.9 on the type-aware model (89.3%). With 30 paired intervals per case size and no correction, about 1.5 would
exclude zero by chance even if Evo 2 added nothing.

: K = 20 (8,230 cases of 21 variants), with 95% intervals. Source: `results.json`, `cases.k20`.

| Measure | Estimate (95% CI) |
|-----------------------------------------------------------------|-----------------------------------:|
| CADD v1.7, in top 5 | 79.6% (76.5 to 82.8) |
| Evo 2 7B, in top 5 | 79.6% (76.3 to 82.4) |
| Type-aware CADD + SpliceAI, in top 5 | 89.3% (87.3 to 91.5) |
| Type-aware CADD + SpliceAI + Evo 2, in top 5 | 90.2% (88.1 to 92.2) |
| Adding Evo 2 to type-aware C + S, all, ranked first | +2.6 points (1.1 to 4.2) |
| Adding Evo 2 to type-aware C + S, all, in top 5 | +0.9 points (0.2 to 1.7) |
| Adding Evo 2 to CADD, coding, ranked first | +2.0 points (-1.2 to 5.2) |
| Adding Evo 2 to CADD, coding, in top 5 | +1.5 points (0.1 to 2.9) |
| Adding Evo 2 to type-aware C + S, noncoding, in top 5 | +0.6 points (0.04 to 1.34) |
| SpliceAI, noncoding, in top 5 | 97.3% (93.7 to 99.3) |
| Type-aware CADD + SpliceAI + Evo 2, noncoding, in top 5 | 93.6% (89.5 to 97.1) |

# Discussion

On one gene with assay labels, Evo 2 alone ranked the loss-of-function variant below CADD v1.7. Added to a
combination of CADD and SpliceAI it raised the top-5 rate by about 3 points. That combination's point estimate was
below CADD alone on coding cases (not tested), and much of the gain may make up for that. Against CADD alone on
coding variants no gain was detected at K = 100. For noncoding variants, no combination's point estimate reached
SpliceAI alone; that comparison was not tested. At K = 20 one coding
comparison against CADD reached nominal significance (+1.5 points in the top 5), which we read as a possible small
effect that this design cannot confirm.

How the scores were combined moved the top-5 rate more than Evo 2 did: from 41.0% to 62.0% without Evo 2,
against gains of 3.1 to 6.0 points from adding it. Evo 2's apparent gain tended to be larger on weaker baselines
(+6.0 points on the average of percentiles against +3.1 on the type-aware model at K = 100; +11.9 against +0.9 at
K = 20). The difference between gains was not tested, and on coding cases at K = 100 the order reversed. A new score
that only helps a weak combination says more about the combination. The rare-disease
benchmark that followed this study used, for that reason, a baseline of all the established scores combined by the
same method [@rewire2026studies; @rewire2026rdv2].

The variant-level comparison of the same scores points the same way [@rewire2026evo2ops]. There, no difference
was detected between Evo 2 and CADD v1.7 overall. CADD v1.7 was better on coding variants: AUROC 0.892 (0.875 to
0.908) against 0.840 (0.818 to 0.861) for Evo 2, a difference of -0.052 (-0.073 to -0.030). On noncoding variants
Evo 2, without training on labels, was not detectably different from SpliceAI: 0.953 (0.935 to 0.968) against 0.963
(0.945 to 0.980), a difference of -0.011 (-0.025 to 0.003). Coming close to a splicing model on splice variants is a
real property of a zero-shot model. It is not a gain for a lab that already runs SpliceAI.

## Limits

- One gene and one assay. Real cases hold rare variants across many genes, after frequency and consequence filters.
  The candidates here come from the same gene's assay, from the same 13 exons as the diagnostic variant.
- Assay labels, not clinical diagnoses. LOF in a cell-based assay is close to, but not the same as, pathogenic.
- Noncoding means almost entirely splice variants here (230 of 239 noncoding LOF variants), so the noncoding result
  is a splice result.
- The logistic models are fitted and tested on the same assay, with position-grouped cross-validation. That keeps
  each position's labels out of its own model but not the gene's.
- REVEL was trained on labelled missense variants that may include BRCA1 [@rewire2026evo2ops]. It scores 8 of the
  1,125 noncoding variants, so it almost never ranks a noncoding LOF variant highly (0.0% in the top 5 at
  K = 100, 0.2% at K = 20).
- Ties count against a score. This is strict for any score with many equal values.
- Intervals resample diagnostic variants by position. They hold fixed the candidate draws, the pool of 3,070
  candidates shared by all cases, and the fitted models, so they understate total uncertainty.
- Logistic scores are out-of-fold predictions from five models, so a case can compare probabilities from different
  fits. Single scores and the average of percentiles are not affected.
- No analysis plan was fixed in advance, and there are 60 paired intervals with no multiplicity correction. All
  comparisons are exploratory. Seven of the 30 paired results at K = 100 change significance at K = 20.
- The run command does not pin package versions. A rerun with the versions in section 2.6 gave identical results.

This benchmark measures ranking signal. It does not validate a diagnostic test.

# Conclusion

In simulated BRCA1 cases, Evo 2 7B alone ranked the loss-of-function variant below CADD v1.7. Added to existing
scores it gave small top-5 gains, which tended to be larger on weaker baselines (not tested). Against CADD alone on
coding variants no gain was detected at 100 candidates per case, and for noncoding variants no combination's point
estimate reached SpliceAI alone (not tested). Whether a new score helps should be tested
against a baseline that is already as good as the existing scores can make it.

# Data and code availability

Code and results: <https://github.com/rewire-bio/amp-2026-studies>, folder `brca1-case-ranking/`.
`case_ranking.py` and `results.json` are unchanged since tag `amp-2026-v1`; this paper, its figures and
`paper/charts.py` were added after tag `amp-2026-v5`. The code is under the MIT licence.

Inputs come from the Evo 2 cloud operations study at commit `f8ae11c` [@rewire2026evo2ops]. Per-variant scores
and labels are not redistributed here, because of their licences:

- Findlay et al. 2018 function scores: free for nonprofit use; commercial use needs a licence
  [@findlay2018accurate].
- CADD v1.7 scores and annotations: free for non-commercial use [@caddteam2026download].
- SpliceAI precomputed scores (Illumina), used as annotated in CADD v1.7 [@illumina2019spliceai].
- REVEL v1.3: free for non-commercial use [@revelteam2026downloads].
- Evo 2 7B weights: Apache-2.0 [@arcinstitute2025evo2].

The folder also holds `liftover_brca1.py` and `fetch_avi.py`, which lift the variants to hg38 and fetch AlphaGenome
AVI scores for them. No AVI results are part of this study.

# References {-}

::: {#refs}
:::
