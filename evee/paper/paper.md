---
title: "How much of a ClinVar AUROC is the consequence class?"
subtitle: "Charts and checks on the EVEE manuscript's published benchmark data"
author: '`Tim Richardson \\ Rewire`{=latex}'
date: 2 November 2026
abstract: |
  Variant effect predictors are usually compared by one area under the ROC curve (AUROC) over ClinVar. EVEE, a
  probe trained on Evo 2's internal states, reports an AUROC of 0.997 on 833,970 ClinVar single-nucleotide variants
  (SNVs), and its authors warn that this figure is largely confounded by consequence type. We redrew two of their
  figures for a talk and checked them against the data they released. We add no new model scores. From the
  published counts, the share of variants labelled pathogenic runs from 0.10% (95% interval 0.09% to 0.11%) for
  synonymous to 98.6% (98.4% to 98.8%) for splice donor variants. A score that knows only the consequence class
  reaches an AUROC of 0.967 (0.967 to 0.968) on the full set. On EVEE's deconfounded set of 158,616 variants it
  falls to 0.737 (0.735 to 0.739). Within each consequence class, AUROCs are lower and the order of methods
  changes. On the deconfounded set EVEE leads CADD v1.7 on intronic variants (0.965 against 0.904) and trails REVEL
  and AlphaMissense on missense variants (0.967 against 0.981 and 0.972). The deconfounded table has no intervals
  and the methods were scored on different subsets, so we give approximate intervals and make no significance claims
  there. On the full set, where EVEE published 95% intervals, EVEE's AUROC is above CADD v1.7's in four of the six
  classes we chart with intervals that do not overlap. Both charts reproduce the published values exactly. The
  practical point for a laboratory: read a predictor's per-consequence results, not its headline AUROC.
---

# Introduction

A variant effect predictor gives each variant a score. Papers usually summarise it with one AUROC over a labelled
set, most often ClinVar [@landrum2018clinvar]. AUROC is the chance that a randomly chosen pathogenic variant scores
higher than a randomly chosen benign one: 0.5 is chance and 1 is perfect.

ClinVar is not a random sample of variants. Whether a variant is labelled pathogenic depends strongly on its
consequence. Nonsense and canonical splice variants are mostly pathogenic; synonymous and intronic variants are
almost all benign. A score that has learned consequence, directly or through features such as conservation, can
reach a high AUROC without telling apart variants of the same type. Lu and colleagues showed that this kind of
heterogeneity inflates reported performance [@lu2026genomic]. Grimm and colleagues described related circularity
in missense benchmarks [@grimm2015evaluation].

A clinical laboratory gets little from that part of the AUROC. Under ACMG/AMP the consequence already enters the
classification on its own terms: a null variant in a loss-of-function gene can meet PVS1, and a computational
prediction adds PP3 or BP4 on top [@richards2015standards]. What the lab needs from a predictor is to separate
variants within a consequence class.

EVEE is a probe trained on the internal states of the Evo 2 DNA model [@brixi2026genome] to predict ClinVar labels
[@pearce2026interpretable]. Its authors report an overall SNV AUROC of 0.997 and say plainly that this figure is
"largely confounded by a model's ability to recover consequence-type priors rather than variant-level
discrimination" [@pearce2026interpretable]. They built a deconfounded subset and compared methods within each
consequence class, and they released the figure data [@goodfire2026evee].

For a talk on genomic foundation models in rare-disease diagnosis we redrew two of their figures in the talk's
colours. This paper documents those charts and the checks we ran on them. It asks three questions:

1. How strongly does consequence predict the ClinVar label in EVEE's sets? We answer with the pathogenic share per
   class and with the AUROC of consequence alone.
2. What do the within-consequence comparisons on the deconfounded set show, and how uncertain are they?
3. Do our charts show the published values?

This is not a new benchmark. Every model score comes from EVEE's released tables. The only values we compute are
derived from their counts and AUROCs, and we mark them as ours.

# Methods

## Data

All inputs are three files in `artifacts/` of the goodfire-ai/evee-manuscript repository at commit `bb3bfa6`
(16 October 2026) [@goodfire2026evee], read from a local clone and not copied into this repository:

- `dataset_characterization.feather`: for each consequence class, the number of variants and the number labelled
  pathogenic, in two sets. "Naive (834K)" is the full SNV set; "Deconf (159K)" is the deconfounded set.
- `deconf_benchmark.feather`: for each method and consequence class on the deconfounded set, the number of variants,
  the number pathogenic, the number the method scored (`n_valid`), AUROC and AUPRC. It has no interval columns.
- `snv_benchmark.feather`: the same for the full SNV set, with a 95% interval for most AUROCs (`auroc_lo`,
  `auroc_hi`; missing for some sparsely scored rows, such as REVEL on synonymous and nonsense variants). The repository does not include the code that computed these intervals, so we do not know the method.

The EVEE repository describes how the sets were made [@goodfire2026evee]. Labels come from ClinVar release
2026-03-21 on GRCh38 with GENCODE v49: pathogenic or likely pathogenic against benign or likely benign, at least
one review star, genes no longer than 100 kb, with multiallelic, multinucleotide and unannotated entries removed.
The deconfounded subset reduces label imbalance within consequence types "while accounting for their prevalence"
(README), with variants selected round-robin over genes; the preprint's Methods give the details
[@pearce2026interpretable]. It is not a scaled-down copy of the full set: missense variants are 15% of the full set
and 64% of the deconfounded set, synonymous variants 35% and 8% (Table 1). The docstring
of the repository's Figure 1c script calls this set "CADD-deconfounded v3". We could not find out from the
repository whether CADD scores were used in selecting it.

## How EVEE scores a variant

Figure 1 summarises the EVEE probe from the repository's protocol [@goodfire2026evee]. Evo 2 7B reads a 65 kb
window around the variant, for the reference and the variant sequence, on both strands. The probe takes the hidden
state of block 27, keeps the 256 positions per strand where reference and variant differ most by cosine divergence,
and forms variant minus reference features. A covariance pooling layer and a linear readout, both trained on
ClinVar SNV labels, give the probability that the variant is pathogenic. Genes were split into five partitions; each
partition was predicted by a probe trained on the other four. The chart label "EVEE probe" is the row named
"Covpool probe" in the tables. The tables also have an "Evo2 probe" row, a mean-pooled probe, which we do not chart.

![How EVEE scores one SNV. Steps 2 and 3 use Evo 2 unchanged; steps 4 and 5 are trained on ClinVar labels with
genes held out. From the evee-manuscript README [@goodfire2026evee].](figures/diagram_probe.pdf){width=85%}

## Methods compared

The results chart shows five of the table's methods: the EVEE probe; Evo 2 zero-shot ("Evo2 loss"), the change in
Evo 2's sequence likelihood when the variant is swapped in [@brixi2026genome]; AlphaMissense
[@cheng2023accurate]; REVEL [@ioannidis2016revel]; and CADD v1.7 [@schubach2024cadd]. It shows six consequence
classes: missense, synonymous, intronic, splice donor, splice acceptor and nonsense. A point is drawn only where the
method scored at least 100 variants in that class. A method that cannot score a variant is left out of that class,
not counted as wrong, so each AUROC describes the variants that method scored.

## The charts

`evee/evee_charts.py` draws both charts with matplotlib 3.11.2, pandas 3.1.0 and pyarrow 26.0.0 under Python
3.11.13, run by uv 0.8.2:

```
uv run --python 3.11 --with matplotlib==3.11.2 --with pandas==3.1.0 \
    --with pyarrow==26.0.0 \
    python evee/evee_charts.py --evee ../evee-manuscript/artifacts
```

The confounding chart (Figure 3) plots the pathogenic share of eight classes in the full set. The results chart
(Figure 4) plots AUROC by class on the deconfounded set. The script writes SVG to `evee/out/`. `paper/make_figures.sh`
reruns it, converts the SVGs to PDF and PNG with rsvg-convert 2.62.3, runs `paper/stats.py`, then renders the
diagrams with mmdc 11.16.0. Figure 2 shows the flow.

![What this folder does with EVEE's released files.](figures/diagram_pipeline.pdf){width=100%}

## Values computed here

`evee/paper/stats.py` (pandas 3.1.0, pyarrow 26.0.0 and numpy 2.4.6, same Python) writes every value below to `evee/paper/stats.json`.
It needs no GPU and runs in seconds on a laptop.

- **Pathogenic share** per class, with a Wilson 95% interval [@wilson1927probable] from the published counts.
- **Consequence alone as a score.** Each variant gets the pathogenic share of its class in that set. The AUROC of
  this score follows exactly from the counts, with ties counted as one half. The 95% interval comes from 2,000
  bootstrap draws (seed 20261102) that resample the pathogenic and the benign variants separately across classes and
  recompute the shares in each draw. The shares are fitted on the same data they are tested on. To measure the
  optimism from this, we also compute a two-fold cross-fitted AUROC: each class's pathogenic and benign counts are
  split at random into halves, the shares are fitted on one half and scored on the other, both ways, and the result
  is averaged over 200 splits (seed 20261103).
- **Pairs within a class and majority-rule accuracy.** From the same counts: the share of pathogenic and benign
  pairs that fall in the same consequence class, and the accuracy of predicting each class's majority label.
- **Approximate intervals for the deconfounded AUROCs**, by the Hanley and McNeil formula
  [@hanley1982meaning] from the number of pathogenic and benign variants each method scored. Where a method scored
  only part of a class, the table does not say how many of the scored variants were pathogenic. We then assume the
  class's pathogenic share and mark the interval as assumed. This affects REVEL, CADD v1.7 on missense and synonymous
  (8 and 1 variants unscored) and the overall rows of AlphaMissense and CADD. For AlphaMissense on missense, the
  scored count equals the number of canonical missense variants, whose split is published, so we use that.
- **Chart check.** The script reads the SVGs back. For the results chart it converts each marker's horizontal
  position to AUROC using the axis ticks at 0.70 and 1.00, and compares it with the table. For the confounding chart
  it compares the printed labels with the counts.

## Statistics

No hypothesis was set in advance; every comparison here is exploratory. On the full SNV set we call a difference
between two methods significant at the 5% level when their published 95% intervals do not overlap. We apply this
rule only between methods that scored the same variants (EVEE, Evo 2 zero-shot and CADD v1.7). For two
independent or positively correlated estimates this rule is stricter than a 5% test. Overlapping intervals do not
show that there is no difference. We make no correction for multiple comparisons. On the deconfounded set we make
no significance claims: the intervals are approximate, and the methods were scored on different, overlapping
subsets, so their estimates are neither paired nor independent.

# Results

## Consequence predicts the label

Table 1 gives the pathogenic share by class in both sets. In the full set, nonsense, splice donor and splice
acceptor variants are 98.2% to 98.6% pathogenic. Synonymous variants are 0.10% pathogenic and intronic variants
0.68%. Intronic and synonymous variants make up 547,106 of the 833,970 SNVs. Figure 3 shows eight of the classes.

| Consequence | Full set: n | Pathogenic, % (95% CI) | Deconfounded: n | Pathogenic, % (95% CI) |
|---|---:|---:|---:|---:|
| Splice donor | 11,855 | 98.6 (98.4 to 98.8) | 1,853 | 92.4 (91.2 to 93.6) |
| Nonsense | 36,050 | 98.3 (98.2 to 98.5) | 6,367 | 91.4 (90.7 to 92.1) |
| Splice acceptor | 9,828 | 98.2 (98.0 to 98.5) | 1,650 | 90.7 (89.2 to 92.0) |
| Initiator codon | 1,363 | 87.4 (85.5 to 89.0) | 715 | 78.2 (75.0 to 81.1) |
| Missense | 125,092 | 29.0 (28.8 to 29.3) | 101,666 | 32.0 (31.7 to 32.3) |
| Other | 3,498 | 7.4 (6.6 to 8.3) | 187 | 54.5 (47.4 to 61.5) |
| Intronic | 254,737 | 0.68 (0.65 to 0.71) | 27,644 | 4.95 (4.70 to 5.21) |
| 5' UTR | 15,232 | 0.54 (0.43 to 0.67) | 1,480 | 5.07 (4.06 to 6.31) |
| 3' UTR | 16,679 | 0.16 (0.11 to 0.24) | 876 | 2.85 (1.94 to 4.18) |
| Non-coding transcript | 67,267 | 0.16 (0.13 to 0.19) | 3,593 | 2.70 (2.22 to 3.28) |
| Synonymous | 292,369 | 0.10 (0.09 to 0.11) | 12,585 | 2.30 (2.05 to 2.57) |
| All | 833,970 | 11.6 | 158,616 | 27.8 |

: Pathogenic share by consequence in EVEE's two ClinVar SNV sets. Counts from `dataset_characterization.feather`;
Wilson intervals computed here. All: 96,793 and 44,057 pathogenic.

A score that knows only the consequence class reaches an AUROC of 0.967 (0.967 to 0.968) on the full set. That is
above Evo 2 zero-shot on the same set (0.940, 0.939 to 0.941). On the deconfounded set, consequence alone gives
0.737 (0.735 to 0.739). Fitting the shares on half the data and scoring the other half gives 0.967 and 0.737 again
(0.96718 and 0.73675 against 0.96720 and 0.73680 in-sample), so the in-sample fit adds less than 0.0001. The
preprint reports a similar fall in consequence-only accuracy, from 96% to 73% (Methods, "Clinical variant
datasets") [@pearce2026interpretable]; its exact rule is not given, and a majority rule per class on the released
counts gives 95.2% and 77.7%. The deconfounding narrows the gaps
between classes but does not remove them: pathogenic shares still run from 2.30% (synonymous) to 92.4% (splice
donor).

![Share of variants labelled pathogenic, by consequence, in EVEE's full set of 833,970 ClinVar SNVs. Drawn from
`dataset_characterization.feather`.](figures/evee_confounding.pdf){width=100%}

## Within each consequence, the deconfounded set

Figure 4 is the results chart and Table 2 gives its values. No method reaches the overall AUROC of 0.997 in any of
these classes. CADD v1.7 ranges from 0.789 (nonsense) to 0.952 (missense) within a class, against 0.963 overall on
this set.

![AUROC within each consequence class on EVEE's deconfounded ClinVar set. Points where a method scored fewer than
100 variants are left out. Drawn from `deconf_benchmark.feather`.](figures/evee_results.pdf){width=100%}

| Class | Method | AUROC | Approximate 95% CI | Variants scored |
|---|---|---:|---:|---:|
| Missense | EVEE probe | 0.967 | 0.964 to 0.971 | 20,679 |
|  | Evo 2 zero-shot | 0.939 | 0.937 to 0.941 | 101,666 |
|  | AlphaMissense | 0.972 | 0.971 to 0.974 | 95,791 |
|  | REVEL | 0.981 | 0.979 to 0.983\* | 29,113 |
|  | CADD v1.7 | 0.952 | 0.951 to 0.954\* | 101,658 |
| Synonymous | EVEE probe | 0.928 | 0.879 to 0.978 | 2,557 |
|  | Evo 2 zero-shot | 0.753 | 0.720 to 0.786 | 12,585 |
|  | CADD v1.7 | 0.933 | 0.913 to 0.953\* | 12,584 |
| Intronic | EVEE probe | 0.965 | 0.951 to 0.980 | 5,771 |
|  | Evo 2 zero-shot | 0.801 | 0.787 to 0.815 | 27,644 |
|  | REVEL | 0.789 | 0.548 to 1.000\* | 106 |
|  | CADD v1.7 | 0.904 | 0.893 to 0.915 | 27,644 |
| Splice donor | EVEE probe | 0.871 | 0.825 to 0.916 | 391 |
|  | Evo 2 zero-shot | 0.756 | 0.722 to 0.791 | 1,853 |
|  | CADD v1.7 | 0.860 | 0.837 to 0.883 | 1,853 |
| Splice acceptor | EVEE probe | 0.907 | 0.871 to 0.943 | 377 |
|  | Evo 2 zero-shot | 0.784 | 0.753 to 0.815 | 1,650 |
|  | CADD v1.7 | 0.865 | 0.843 to 0.887 | 1,650 |
| Nonsense | EVEE probe | 0.823 | 0.793 to 0.853 | 1,320 |
|  | Evo 2 zero-shot | 0.784 | 0.768 to 0.800 | 6,367 |
|  | CADD v1.7 | 0.789 | 0.773 to 0.805 | 6,367 |

: AUROC on the deconfounded set, from `deconf_benchmark.feather`, with approximate 95% intervals computed here
(Hanley and McNeil). \* The pathogenic count among scored variants is
not published; the interval assumes the class's share. Methods with fewer than 100 variants scored in a class are left out, as on the chart.

What the table shows, with no significance claims (see Statistics):

- **Intronic.** EVEE 0.965 against CADD v1.7 0.904. The approximate intervals are well apart. This is the largest
  gap to CADD in the chart.
- **Missense.** EVEE 0.967 is below REVEL (0.981) and AlphaMissense (0.972) and above CADD v1.7 (0.952). The
  approximate intervals are apart for REVEL and CADD v1.7. For AlphaMissense they almost touch (EVEE up to 0.9706,
  AlphaMissense from 0.9710). EVEE's missense AUROC is on 20,679 variants,
  AlphaMissense's on 95,791 and REVEL's on 29,113.
- **Splice donor, splice acceptor, nonsense, synonymous.** EVEE is above CADD v1.7 in the first three and just
  below it on synonymous variants (0.928 against 0.933). Every one of these pairs of intervals overlaps. EVEE's
  subsets here are small: 391, 377, 1,320 and 2,557 variants, with only 50 pathogenic synonymous variants.
- **EVEE against Evo 2 zero-shot.** The probe is higher in all six classes. The intervals are apart in five; for
  nonsense they overlap (0.793 to 0.853 against 0.768 to 0.800).
- **REVEL on intronic variants** rests on 106 scored variants, about 5 of them pathogenic if the class share holds.
  Its interval runs from 0.548 to 1.000. The point on the chart carries almost no information.

The methods were not scored on the same variants. EVEE was scored on 32,437 of the 158,616 deconfounded variants,
20.5% of them; CADD v1.7 on 158,607 and Evo 2 zero-shot on all 158,616. On the full set, by contrast, EVEE has a
score for all 833,970 variants. The table's Overall row gives 833,970 scored variants, but its "Other" row for the
EVEE probe has 390 variants against 3,498 for every other method, so its class rows sum to 830,862. The six classes
we report are complete. The repository does not explain the difference. A fraction of about one in five
would fit a single held-out gene partition, but we could not confirm this.

## Within each consequence, the full set

The full-set table has published 95% intervals and EVEE scores every variant in the six classes we report, so it
allows comparisons on nearly the same variants for EVEE, Evo 2 zero-shot and CADD v1.7; AlphaMissense and REVEL
scored subsets. Table 3 gives the same five methods and six classes.

| Class | Method | AUROC | Published 95% CI | Variants scored |
|---|---|---:|---:|---:|
| All SNVs | EVEE probe | 0.9968 | 0.9967 to 0.9970 | 833,970 |
|  | Evo 2 zero-shot | 0.9401 | 0.9393 to 0.9410 | 833,970 |
|  | CADD v1.7 | 0.9914 | 0.9911 to 0.9917 | 832,690 |
| Missense | EVEE probe | 0.971 | 0.970 to 0.972 | 125,092 |
|  | Evo 2 zero-shot | 0.926 | 0.924 to 0.927 | 125,092 |
|  | AlphaMissense | 0.968 | 0.967 to 0.969 | 117,653 |
|  | REVEL | 0.978 | 0.977 to 0.980 | 36,783 |
|  | CADD v1.7 | 0.946 | 0.945 to 0.947 | 125,082 |
| Synonymous | EVEE probe | 0.961 | 0.948 to 0.973 | 292,369 |
|  | Evo 2 zero-shot | 0.670 | 0.643 to 0.700 | 292,369 |
|  | CADD v1.7 | 0.929 | 0.909 to 0.950 | 292,368 |
| Intronic | EVEE probe | 0.984 | 0.980 to 0.987 | 254,737 |
|  | Evo 2 zero-shot | 0.801 | 0.787 to 0.814 | 254,737 |
|  | REVEL | 0.771 | 0.604 to 0.908 | 752 |
|  | CADD v1.7 | 0.910 | 0.900 to 0.918 | 254,736 |
| Splice donor | EVEE probe | 0.920 | 0.894 to 0.944 | 11,855 |
|  | Evo 2 zero-shot | 0.762 | 0.714 to 0.811 | 11,855 |
|  | CADD v1.7 | 0.857 | 0.819 to 0.893 | 11,855 |
| Splice acceptor | EVEE probe | 0.928 | 0.902 to 0.951 | 9,828 |
|  | Evo 2 zero-shot | 0.800 | 0.759 to 0.840 | 9,828 |
|  | CADD v1.7 | 0.882 | 0.841 to 0.916 | 9,828 |
| Nonsense | EVEE probe | 0.900 | 0.885 to 0.914 | 36,050 |
|  | Evo 2 zero-shot | 0.809 | 0.788 to 0.828 | 36,050 |
|  | CADD v1.7 | 0.783 | 0.760 to 0.801 | 36,050 |

: AUROC on EVEE's full SNV set with the published 95% intervals, from `snv_benchmark.feather`. EVEE and Evo 2
zero-shot scored all variants in each class; CADD v1.7 all but 1,280 of 833,970, nearly all of them in the
"Other" class. Rows with fewer than 100 scored variants (AlphaMissense on splice acceptor, 29) and rows without a
published interval (REVEL on synonymous and nonsense variants) are left out.

By the rule in Statistics, all exploratory:

- **EVEE against CADD v1.7.** EVEE is higher overall (0.997 against 0.991) and within missense, intronic, splice
  donor and nonsense variants, with intervals that do not overlap. For splice donor the gap between the intervals is
  0.001. On synonymous and splice acceptor variants the intervals overlap: not shown.
- **EVEE against Evo 2 zero-shot.** EVEE is higher overall and in all six classes, with intervals that do not
  overlap.
- **Missense against the protein tools.** EVEE has 0.971, AlphaMissense 0.968 and REVEL 0.978. AlphaMissense and
  REVEL scored 117,653 and 36,783 of the 125,092 missense variants, so these AUROCs are on different variants and we
  do not apply the overlap rule to them. REVEL is higher on the missense variants it scores.

The overall AUROC compresses large differences. EVEE and CADD v1.7 differ by 0.005 overall, but by 0.025 on
missense, 0.074 on intronic and 0.117 on nonsense variants. CADD v1.7 ranges from 0.783 to 0.946 across these six
classes against 0.991 overall.

EVEE's within-class AUROCs are higher on the full set than on the deconfounded set: intronic 0.984 against 0.965,
splice donor 0.920 against 0.871, nonsense 0.900 against 0.823. CADD v1.7 changes less between the two (intronic
0.910 against 0.904, nonsense 0.783 against 0.789). The two sets differ in composition and EVEE's deconfounded
subset is smaller, so we do not interpret this further.

## The charts match the published values

All 21 markers on the results chart decode to the table's AUROC; the largest difference is $2 \times 10^{-10}$
(`stats.json`, `max_results_diff`). All 8 labels on the confounding chart match the counts to the one decimal place
shown.

# Discussion

Consequence alone reaches an AUROC of 0.967 on EVEE's full ClinVar SNV set. In the full set 95% of
pathogenic-benign pairs are in different consequence classes (within-class pairs are 5.3%), so the overall AUROC
mostly measures how a score ranks classes. So a headline AUROC of 0.99 or more says little about how well a tool separates
two missense variants, or two intronic variants. The overall gap between EVEE and CADD v1.7 is 0.005; within a
class it is 0.025 to 0.117 on the full set. Within-class comparisons, which the EVEE authors report themselves,
are the numbers a laboratory should read.

The within-class results favour EVEE most outside missense variants. On intronic variants EVEE's lead over CADD
v1.7 is large on both sets. On missense variants REVEL is higher on both sets, on the variants it scores (36,783 of
125,092 on the full set), and AlphaMissense is within 0.005 of it on both. The classes outside missense are where a DNA model could add something that protein tools cannot: AlphaMissense has no intronic,
synonymous or splice rows at all on the deconfounded set, and REVEL scores only 106 of 27,644 intronic variants.

The probe also beats Evo 2's own zero-shot score in every class. The EVEE authors take this as evidence that
Evo 2's internal states carry more than its likelihood does [@pearce2026interpretable]. Part of the gain may come
from the labels: the probe is trained on ClinVar and tested on ClinVar, in held-out genes, while CADD v1.7 and the
zero-shot score never saw ClinVar labels. The authors note that ClinVar benchmarks favour models trained on ClinVar
[@pearce2026interpretable].

## Limits

- **No new scores.** We did not run any model. Every AUROC is EVEE's. Our own work is the charts, the check that
  they match, and values derived from published counts.
- **Inputs not in this repository.** The three feather files are read from a clone of goodfire-ai/evee-manuscript
  at commit `bb3bfa6`. If that repository changes, the charts can change with it.
- **No published intervals on the deconfounded set.** Our Hanley and McNeil intervals are approximate. The formula
  assumes a particular score distribution, and for REVEL and part of CADD the pathogenic count among scored
  variants is assumed. We make no significance claims on this set.
- **Different subsets.** On the deconfounded set EVEE was scored on 20.5% of the variants and the other methods on
  more. A tool that cannot score a variant is left out, not counted as wrong. Two AUROCs in the same row of Table 2
  are therefore on different variants. We do not know why EVEE has a full-set score for every variant but a
  deconfounded score for only a fifth.
- **Unknown interval method on the full set.** The published intervals come without the code that made them.
- **The deconfounded set's construction.** It is described in the preprint, not in code we could run, and its
  script names it "CADD-deconfounded". If CADD scores were used in selecting it, CADD's AUROCs on it could be biased.
- **Consequence-only AUROC is in-sample.** The shares are fitted on the set they score. The two-fold cross-fitted
  value is lower by less than 0.0001 on both sets.
- **Exploratory comparisons.** None were planned in advance, and there is no correction for multiple comparisons.
- **ClinVar labels.** The labels are themselves partly based on computational evidence and on consequence, which
  can favour tools that use the same evidence [@grimm2015evaluation].

# Conclusion

In EVEE's ClinVar SNV set the consequence class alone gives an AUROC of 0.967. Headline AUROCs near 0.99 mostly
measure that. Within each class, AUROCs are lower and the order of methods changes: EVEE leads clearly on intronic
variants, REVEL is higher on the missense variants it can score (36,783 of 125,092 on the full set), and on the
deconfounded set EVEE's gaps to CADD v1.7 in the synonymous, splice and nonsense classes are within the approximate
intervals. Our two charts show EVEE's published values exactly. A lab choosing a predictor should look at its
results within the consequence classes that matter for its cases, on variants the tool can actually score.

# Data and code availability

- **Code.** `evee/evee_charts.py` (the charts), `evee/paper/stats.py` (the values computed here, written to
  `evee/paper/stats.json`) and `evee/paper/make_figures.sh` are in the repository rewire-bio/amp-2026-studies
  [@rewire2026studies], under the MIT licence. The chart script has been in every release since tag
  `amp-2026-v4`. This paper, `stats.py`, `make_figures.sh` and the shared build tooling in
  `tools/paper/` are not yet in a tagged release; the latest tag is `amp-2026-v5`.
- **Rebuild.** With goodfire-ai/evee-manuscript cloned next to the repository at commit `bb3bfa6`, run
  `bash evee/paper/make_figures.sh`, then `bash tools/paper/build.sh evee/paper`.
- **Inputs.** `dataset_characterization.feather`, `deconf_benchmark.feather` and `snv_benchmark.feather` from
  goodfire-ai/evee-manuscript at commit `bb3bfa6` [@goodfire2026evee]. The repository has no licence file at that
  commit, so we do not copy the files; the counts and AUROCs we use are in
`evee/paper/stats.json`. The per-variant EVEE scores on Zenodo (10.5281/zenodo.19701997) are
  under CC BY 4.0, but we did not use them here. ClinVar, the source of the labels, is in the public domain
  [@landrum2018clinvar].

