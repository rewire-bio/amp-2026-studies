---
title: "Evo 2 variant calls change between GPU models"
subtitle: "BRCA1 scores on H100, L40S and L40 GPUs with everything else held fixed"
author: '`Tim Richardson \\ Rewire`{=latex}'
date: 2 November 2026
abstract: |
  A lab that runs a DNA language model itself has to choose a GPU. We asked whether that choice changes the
  model's output. The same container image, Evo 2 7B weights and inputs scored 3,893 BRCA1 single-nucleotide
  variants with saturation genome editing labels on two NVIDIA H100 GPUs at two cloud providers and on an NVIDIA
  L40S. The two H100s gave bit-identical scores, and every platform reproduced its own scores exactly on fresh
  machines. The L40S did not: per-variant scores differed by a median of $1.17 \times 10^{-4}$ (Spearman 0.983),
  while AUROC changed by only 0.0021 (95% interval 0.0002 to 0.0042). With decision cut-offs fitted on the H100
  scores and applied unchanged to the L40S scores, 72 of 3,893 calls changed (1.8%, exact 95% interval 1.4% to
  2.3%) at the Youden cut-off, and 152 (3.9%, 3.3% to 4.6%) and 189 (4.9%, 4.2% to 5.6%) at cut-offs set for 90% and
  95% sensitivity. Same-GPU controls changed no calls. In an exploratory follow-up on 487 variants, an L40 and an
  L40S (same architecture) disagreed as much as an H100 and an L40S; turning off FP8 alone narrowed the gap to 0.70
  of its default size; and computing in FP32 brought the two GPUs to a median difference of
  $4.8 \times 10^{-7}$ with 0, 0 and 1 changed calls, at 3.4 to 5.8 times the GPU time. A validated Evo 2
  pipeline should fix the GPU product and validate on it, or compute in FP32.
---

# Introduction

Genomic foundation models such as Evo 2 [@brixi2026genome] score a variant by running a neural network, not by
looking it up in a table. A lab that wants an Evo 2 score for a variant that no precomputed table covers has to run
the model on its own GPU (graphics processing unit). Clinical tests are validated as a fixed method: the same input,
software and settings should give the same output. For a model on a GPU, the hardware is one of the settings.

Floating-point arithmetic is not associative. The order in which partial sums are added changes the last bits of a
result, and GPU libraries choose kernels and summation orders for each device. NVIDIA's cuBLAS library promises
bit-wise identical results only on GPUs with the same architecture and the same number of multiprocessors, with one
toolkit version [@nvidia2026cublas]. cuBLAS is only one of the libraries a model calls. Evo 2 7B runs mostly in
bfloat16 and, in Arc Institute's default configuration, uses FP8 input projections and flash attention
[@arcinstitute2026evo2; @micikevicius2022fp8; @dao2022flashattention]. FP8 reduces precision further, and flash
attention changes the order in which partial sums are accumulated. So scores may differ between GPU models. What a lab needs to know is by how much,
and whether the difference changes a call.

This paper reports a study that measured this for Evo 2 7B on BRCA1 [@rewire2026evo2ops]. It asks four questions:

1. Do two GPUs of the same model, at two cloud providers, give the same scores? Does each platform reproduce its
   own scores on fresh machines?
2. How much do scores differ between an H100 and an L40S, and does discrimination against assay labels change?
3. At a decision cut-off fixed on the H100, how many variants change call on the L40S?
4. Exploratory: what causes the difference, and what removes it? FP8, full 32-bit precision (FP32), GPU
   architecture and batch size were tested.

The study was written for a talk to clinical scientists. The charts in this paper show one variant whose call
changes, and all variants at once (Figures 3 and 4). Figure 1 shows the design.

![Study design. Everything except the GPU and, in the follow-up, the numerical configuration was held fixed.](figures/diagram_design.pdf){width=62%}

# Methods

## Source of the data and code

This paper does not run new GPU jobs. All scores, receipts and analysis outputs come from the committed evidence in
the repository `rewire-bio/evo2-cloud-operations` at commit `f8ae11c9a1c4` [@rewire2026evo2ops]: the recorded run
`b75c67aa` (9 October 2026), its independent reproduction `30e7ae49` (same day, fresh machines), the follow-up of
11 October 2026 (protocol amendment 5) and the classification-flip analysis of 13 October 2026 (amendment 6). That
repository holds the protocol, its seven amendments, the launchers and the run logs, including failed attempts. The
charts here are drawn from those files by `gpu-drift/flip_charts.py` and `gpu-drift/paper/make_figures.py` in
`rewire-bio/amp-2026-studies`.

## Data

The variants are the 3,893 BRCA1 single-nucleotide variants (SNVs) of Findlay et al., measured by saturation
genome editing in a human cell line [@findlay2018accurate]: 2,821 functional, 249 intermediate and 823 loss of
function (LOF). 2,768 are coding and 1,125 noncoding. Coding means missense, nonsense or synonymous; noncoding
means canonical splice, splice region, intronic or 5' UTR. The variant table, GRCh37 chromosome 17 and Arc's four
forward-test sequences were fetched from Arc's Evo 2 repository at commit `53f1959` [@arcinstitute2026evo2] and
checked by SHA-256. No patient data were used.

## Model and software

All runs used one container image, `ghcr.io/rewire-bio/evo2-cloud-operations` at digest
`sha256:22970720ff8bc05a00460a6066a70ed71951001aa52965141591a0ca6dcc968b`, built on NGC PyTorch 25.04, with
Python 3.12.3, PyTorch 2.7.0a0+79aa17489c.nv25.04, CUDA 12.9, cuDNN 9.9.0, evo2 0.6.0, vtx 1.1.0, flash-attn 2.7.3
and the Transformer Engine build shipped in that image. The `evo2_7b` weights were pinned to Hugging Face revision
`bda0089f` and checked by SHA-256 (`c6664592`, first 8 characters) before loading. The default configuration is
Arc's: parameters mostly in bfloat16 (some in FP32), FP8 input projections through Transformer Engine, flash
attention, and TF32 matrix arithmetic not allowed. Study code was copied to each machine from revisions `e93491ac`
(recorded run), `64de22d0` (reproduction) and `7ede2bad` (follow-up) of `evo2-cloud-operations`.

## Scoring

Each SNV was scored as in Arc's BRCA1 notebook [@arcinstitute2026evo2]. An 8,192 bp window centred on the variant
is scored with and without the variant. The score used here is the reference window's mean log-likelihood minus the
variant window's, so higher means more damaging. The 3,893 variants share 1,326 reference windows, so 5,219 windows
were scored, one at a time (batch size 1). Before scoring, each machine ran a reimplementation of Arc's forward
test (workload W0): the mean loss on Arc's four test sequences must be within 0.001 of 0.3477.

Each machine also ran two side workloads. W2 embeds Arc's exon-classifier sample positions with `evo2_7b_base`
(Hugging Face revision `074097e9`) and scores them with the pinned exon classifier
(`schmojo/evo2-exon-classifier`, revision `3ce7ccd4`). W4 reloads `evo2_7b` in a new process and rescores the first
200 BRCA1 variants.

## Compute

Table 1 lists the machines. Each was created fresh with a 3-hour limit, ran the pinned image, copied results back
and was deleted. Runpod pods were on-demand Secure Cloud pods; the Google Cloud machine was Spot.

| Run | Platform | GPU | Driver | Location | USD/h |
|-------------|---------------------|------|-----------|-------------|------:|
| Recorded | P2, Google Cloud | H100 | 580.178.04 | us-central1-a | 6.24 |
| Recorded | P3, Runpod | H100 | 580.126.09 | US-MO-1 | 3.99 |
| Recorded | P5, Runpod | L40S | 580.159.03 | US-MO-1 | 1.09 |
| Reproduction | P2, Google Cloud | H100 | 580.178.04 | us-central1-c | 6.24 |
| Reproduction | P3, Runpod | H100 | 580.126.09 | AP-IN-1 | 3.99 |
| Reproduction | P5, Runpod | L40S | 580.159.03 | US-MO-1 | 1.09 |
| Follow-up | FH, Runpod | H100 | 580.126.09 | US-MO-1 | 3.99 |
| Follow-up | FL, Runpod | L40S | 595.91.07 | US-TX-4 | 1.09 |
| Follow-up | FA, Runpod | L40 | 595.91.07 | EUR-IS-2 | 0.82 |

: Machines used. H100: NVIDIA H100 80GB HBM3, Hopper architecture, compute capability 9.0. L40S and L40: 48 GB,
Ada Lovelace architecture, compute capability 8.9 [@nvidia2026h100; @nvidia2026l40s; @nvidia2026l40]. P2 is a
Google Cloud `a3-highgpu-1g` Spot machine. Prices are list prices recorded at launch. Drivers and locations are from
each machine's receipt.

The recorded run cost USD 4.78 on P2 (2,734 billed seconds), USD 2.82 on P3 and USD 3.78 on P5 over two attempts
(the first failed when an SSH session dropped). The reproduction cost USD 4.62 on P2, USD 2.74 on P3 and USD 2.11
on P5. The follow-up cost USD 4.98. The protocol also planned two AWS platforms (P1 H100, P4 L40S); amendment 2
deferred them because no AWS account was available.

## Follow-up configurations (exploratory)

The follow-up scored every 8th BRCA1 SNV in table order (487 variants, 104 LOF) on one H100, one L40S and one L40,
each configuration in a fresh process. Configurations: `default`; `fp8off` (FP8 input projections off, otherwise
default); `fp32` (all parameters in FP32, TF32 off, PyTorch attention in place of flash attention); and on the H100
only, `default` with batch size 8. These configurations were set by patching Evo 2 at load time
(`followup/gpu_followup.py`): `fp8off` sets `evo2.models.HAS_TE = False`, so input projections run in bfloat16
through Transformer Engine; `fp32` also overrides the block dtypes to float32, disables flash attention and casts
the model to float32. TF32 was disabled for matmul and cuDNN in every configuration. The `default` batch-1 runs on the H100 and L40S were controls: they had to
reproduce the recorded run's scores for the subset exactly. Amendment 5 fixed the analysis before the run and
labelled the follow-up exploratory.

## Agreement statistics

For each pair of score sets: the median and maximum absolute per-variant difference; Spearman correlation; and the
number of variants whose rank moves by more than 1% of the list (39 or more places for 3,893 variants; 5 or more
for the 487-variant subset). Discrimination is the AUROC for LOF against functional or intermediate variants, with
95% percentile intervals from 2,000 bootstrap resamples of positions (all SNVs at a position kept together; 1,326
positions), seed 7. The paired L40S minus H100 difference used a separate set of 2,000 paired resamples (seed 11),
computed in `evidence/checks/verify_results.py`. This difference was computed in the independent checks after the methods review, not in the original protocol.

## Changes of call at a fixed cut-off

Amendment 6 was added on 13 October 2026, after the scores and agreement statistics were known, in response to an
audience question. It fixed the cut-off rules and counts before any flip count was computed. A variant is called LOF if its score is at
or above a cut-off $t$. Candidate cut-offs are midpoints between adjacent distinct reference scores. Three cut-offs
were fitted on the Google Cloud H100 scores (P2, recorded run): the Youden cut-off, which maximises sensitivity plus
specificity minus one [@youden1950index] (ties to the lowest), and the highest cut-offs with sensitivity of at least
0.90 and 0.95. Arc's notebook reports only AUROC and defines no cut-off. Each cut-off was applied unchanged to the
other platform's scores, and a variant whose call differs is a flip.

For each cut-off we report flips overall and by stratum (LOF and not LOF; coding and noncoding), the direction
(LOF on the H100 only, or on the other GPU only), and an exact Clopper-Pearson 95% interval on the proportion
[@clopper1934use]. The distance of each flipped variant from the cut-off is expressed in units of the pair's median
absolute difference. Controls that had to give 0 flips: P2 against P3, and each platform's recorded run against its
reproduction. On the follow-up subset, cut-offs were refitted on the reference member's own subset scores (H100, or
L40 for the L40 pairs), because changing precision shifts every score by about $10^{-4}$.

## Which results are tests

Amendment 6 specified counts with exact intervals, not hypothesis tests. Scoring is deterministic: every platform
reproduced its own scores bit for bit, so any nonzero flip count already shows the rate is not zero. The exact
intervals give the range of flip rates consistent with the counts for comparable variant sets; the comparison that
matters is with the same-GPU controls, which gave exactly 0. The AUROC difference is significant at the 5% level
if its paired bootstrap interval excludes zero. Comparisons between strata, the direction of flips and every
follow-up comparison were not tested and are descriptive. The follow-up is exploratory throughout, and its subset
cut-offs rest on 104 LOF variants.

## Figures

`make_figures.sh` in the paper folder draws the charts with matplotlib from the files above and the diagrams with
Mermaid. It needs a checkout of `evo2-cloud-operations` next to this repository.

![How one call flips. The score is computed on each GPU; the cut-off is fitted on the H100 and applied unchanged to
the L40S.](figures/diagram_flip.pdf){width=48%}

# Results

## Same GPU model: identical scores

Both H100s gave the same score for all 3,893 variants (maximum absolute difference 0, Spearman 1.000), and the same
probabilities in the exon classifier workload (W2). Each platform reproduced its own scores exactly in a new process
on the same machine (W4, 200 variants) and in the reproduction run on fresh machines (all 3,893). The reproduction's
Runpod H100 ran in India (AP-IN-1) rather than Missouri. In the follow-up, the L40S control ran on a different host,
in Texas, with a different driver (595.91.07 against 580.159.03), and still reproduced the recorded L40S scores
for the subset exactly. These are exact equalities, so no test applies.

## H100 against L40S: scores differ, AUROC barely moves

All three platforms passed the forward test with mean loss 0.34765625. The test has bfloat16 resolution: token
accuracy differed on the L40S (0.8626 against 0.8636), but the loss did not. So the check that the installation
works does not detect differences of the size that follow.

The L40S scores differed from the H100 scores by a median of $1.17 \times 10^{-4}$ and at most
$1.31 \times 10^{-3}$. Spearman correlation was 0.983, and 2,627 of 3,893 variants moved by 39 or more rank
places (more than 1% of the list). Among the 823 LOF variants the correlation was 0.998. Of the 50 highest-scoring
variants overall, 48 were shared (98 of 100, 195 of 200).

AUROC was 0.875 (95% interval 0.859 to 0.891) on the H100 and 0.877 on the L40S. The paired difference, L40S minus
H100, was 0.0021 (0.0002 to 0.0042). It is significant at the 5% level and too small to matter for discrimination.

## Calls that change at a fixed cut-off

The Youden cut-off fitted on the H100 was 0.001698, with sensitivity 0.793 and specificity 0.782 on the H100 scores.
The 90% and 95% sensitivity cut-offs were 0.000997 (specificity 0.643) and 0.000698 (specificity 0.539).

Figure 3 follows one variant, chr17:41276094:C>A, a missense variant. On the H100 it scores 0.002101, above the
Youden cut-off, so it is called LOF. On the L40S it scores 0.001625, below the cut-off, so it is called functional.
The score moved by $4.76 \times 10^{-4}$. The assay classes this variant as functional, so here the L40S call agrees with the assay; across all
variants, flips run both ways. It was chosen as the example because, of the variants called LOF on the
H100 and functional on the L40S, its score moved most.

![One variant on both GPUs, around the Youden cut-off (dashed). Blue: H100 score. Orange: L40S score.](figures/gpu_flip_step2.pdf){width=100%}

![All 3,893 variants, H100 score against L40S score, zoomed to the region around the Youden cut-off (dashed lines).
Orange: the 72 variants whose call differs.](figures/gpu_flip_scatter.pdf){width=62%}

Table 2 gives the flips. Flips occurred at every cut-off and in every stratum; the exact intervals give the range
of flip rates consistent with these counts for comparable variant sets. The four controls gave 0 flips at every cut-off, and the reproduction's
L40S gave the same flips as the recorded one.

| Cut-off | All | LOF | Not LOF | Coding | Noncoding | LOF on one GPU only (H100 / L40S) |
|---------|----------|----------|----------|----------|----------|-----------:|
| Youden | 72, 1.8% (1.4\ to\ 2.3) | 13, 1.6% (0.8\ to\ 2.7) | 59, 1.9% (1.5\ to\ 2.5) | 63, 2.3% (1.8\ to\ 2.9) | 9, 0.8% (0.4\ to\ 1.5) | 41 / 31 |
| 90% sensitivity | 152, 3.9% (3.3\ to\ 4.6) | 18, 2.2% (1.3\ to\ 3.4) | 134, 4.4% (3.7\ to\ 5.1) | 125, 4.5% (3.8\ to\ 5.4) | 27, 2.4% (1.6\ to\ 3.5) | 73 / 79 |
| 95% sensitivity | 189, 4.9% (4.2\ to\ 5.6) | 15, 1.8% (1.0\ to\ 3.0) | 174, 5.7% (4.9\ to\ 6.5) | 152, 5.5% (4.7\ to\ 6.4) | 37, 3.3% (2.3\ to\ 4.5) | 102 / 87 |

: Calls that change between the H100 and the L40S at cut-offs fitted on the H100. Count, percentage and exact 95%
interval in percent. Group sizes: all 3,893; LOF 823; not LOF 3,070; coding 2,768; noncoding 1,125.

More calls flip at the lower cut-offs because more variants lie near them: 150, 235 and 334 H100 scores lie within
one median GPU difference of the Youden, 90% and 95% sensitivity cut-offs. Flipped variants sat close to the
cut-off, at a median of 0.68, 0.67 and 0.66 median differences, and at most 5.4, 5.4 and 3.9.

Descriptive comparisons, not tested: flips ran in both directions in similar numbers (41 and 31 at Youden), so the
L40S was not systematically more or less likely to call LOF. Coding variants flipped more often than noncoding ones
at the Youden cut-off (2.3% against 0.8%; the intervals do not overlap). At the 90% and 95% cut-offs, not-LOF
variants flipped more often than LOF variants. These patterns follow where each group's scores sit relative to the
cut-off and were not pre-specified as questions.

## Follow-up: cause and control (exploratory)

The H100 and L40S default controls reproduced the recorded run's subset scores exactly, and every configuration
passed the forward test. Figure 5 and Table 3 summarise agreement; Figure 6 shows flips.

| Comparison (487 variants) | Median abs. difference | Max abs. difference | Spearman | Rank moves > 1% |
|--------------------------|-----------:|-----------:|--------:|------------:|
| H100 vs L40S, default | $1.19 \times 10^{-4}$ | $5.81 \times 10^{-4}$ | 0.983 | 327 |
| H100 vs L40S, FP8 off | $8.35 \times 10^{-5}$ | $3.52 \times 10^{-4}$ | 0.990 | 295 |
| H100 vs L40S, FP32 | $4.77 \times 10^{-7}$ | $2.92 \times 10^{-6}$ | 1.000 | 0 |
| L40 vs L40S, default | $1.13 \times 10^{-4}$ | $6.96 \times 10^{-4}$ | 0.980 | 320 |
| L40 vs L40S, FP8 off | $7.59 \times 10^{-5}$ | $4.15 \times 10^{-4}$ | 0.990 | 298 |
| H100, batch 8 vs batch 1 | 0 | 0 | 1.000 | 0 |
| H100, FP8 off vs default | $1.43 \times 10^{-4}$ | $1.33 \times 10^{-3}$ | 0.976 | 349 |
| H100, FP32 vs default | $1.19 \times 10^{-4}$ | $1.38 \times 10^{-3}$ | 0.981 | 348 |
| L40S, FP8 off vs default | $1.37 \times 10^{-4}$ | $1.66 \times 10^{-3}$ | 0.977 | 349 |
| L40S, FP32 vs default | $1.16 \times 10^{-4}$ | $1.53 \times 10^{-3}$ | 0.981 | 331 |

: Agreement in the follow-up subset. Spearman for FP32 is 0.999997. Rank moves > 1%: 5 or more places.

FP8 is a cause but not the main one. Turning off FP8 input projections reduced the median H100 against L40S
difference to 0.70 of its default value. Amendment 5 set the bar for "main cause" at below 0.10.

FP32 removes nearly all of the difference. With all parameters in FP32, TF32 off and PyTorch attention, the median
difference fell to $4.77 \times 10^{-7}$, about 250 times smaller, and no variant moved by more than 1% of the list
(Spearman 0.999997). The price was speed: the H100
scored 0.50 windows per second against 2.92 by default (5.8 times slower) and the L40S 0.27 against 0.92 (3.4 times
slower). Peak memory was 29.6 GiB on both, so FP32 fits on a 48 GB GPU for 8 kb windows.

Same architecture is not enough. The L40 and L40S, both Ada Lovelace with compute capability 8.9, disagreed by a
median of $1.13 \times 10^{-4}$, about as much as the H100 and L40S ($1.19 \times 10^{-4}$ on the same subset). Batch size 8 gave bit-identical scores to batch
size 1 on the H100.

Within one GPU, changing precision moved scores as much as changing GPU (Table 3, last four rows). Default scores
therefore carry rounding error of about $10^{-4}$ whichever GPU computes them. Subset AUROC was 0.889 to 0.893 in
every configuration; 104 LOF variants are too few to separate them.

![Absolute per-variant score differences in the follow-up subset (log scale). Dot: median; bar end: maximum. Orange:
between GPUs in reduced precision. Blue: between GPUs in FP32, and batch 8 against batch 1. Grey: precision changes
on one GPU.](figures/differences.pdf){width=90%}

At cut-offs refitted on the subset, the default H100 against L40S pair flipped 8, 13 and 22 of 487 calls at the
Youden, 90% and 95% sensitivity cut-offs: 1.6% (0.7 to 3.2), 2.7% (1.4 to 4.5) and 4.5% (2.9 to 6.8). With the
full-set H100 cut-offs instead, it flipped 8, 18 and 24 (1.6%, 3.7% and 4.9%), consistent with the full-set rates.
In FP32 it flipped 0 (0.0 to 0.8), 0 (0.0 to 0.8) and 1 (0.0 to 1.1). The one FP32 flip sat 1.4 FP32 median
differences, about $6.9 \times 10^{-7}$, from its cut-off. FP8 off gave 3, 14 and 11 flips; the L40 against the
L40S gave 8, 14 and 25 by default and 7, 13 and 12 with FP8 off. At the 90% and 95% cut-offs the FP32 and default
intervals do not overlap; no paired test was specified, so this is descriptive.

![Calls that change, with exact 95% intervals, at each cut-off. Subset rows use cut-offs refitted on the subset's
reference GPU. Orange: pairs computed in reduced precision. Blue: FP32 and the same-model control.](figures/flips_by_threshold.pdf){width=100%}

# Discussion

The GPU model is part of the method. With the image, weights and inputs fixed, two H100s at two providers, on
different drivers, gave identical Evo 2 7B scores, and so did the reproduction's Runpod H100 on a host in India.
Every platform reproduced itself. Changing the GPU model changed scores by about $10^{-4}$. That is small against the spread of the scores
and leaves AUROC nearly unchanged. It is not small against a cut-off: 1.8% to 4.9% of BRCA1 calls changed,
depending on where the cut-off sits. A call hides that the variant was borderline.

The cause is in the arithmetic, not the installation. All GPUs passed Arc's forward test. The L40 and L40S have
the same architecture and the same core counts (18,176 CUDA cores and 142 RT cores, one per multiprocessor, so 142
multiprocessors) [@nvidia2026l40datasheet; @nvidia2026l40s], and in the follow-up they ran one image on one driver
(595.91.07). They therefore meet the conditions under which NVIDIA promises identical cuBLAS results
[@nvidia2026cublas], yet they disagreed. The difference probably comes from other kernels whose choice can depend on
the device, such as those in Transformer Engine, flash attention or Triton. The source repository's follow-up notes
attribute it to differing multiprocessor counts, but the data sheets give the same counts. We did not trace individual
kernels, so this is the likely explanation, not a demonstrated one. FP8 explains part of the difference. FP32
removes nearly all of it.

For a lab, three controls follow. Fix the GPU product, not just the architecture, and validate on it; revalidate
when it changes. If scores must depend as little as possible on the hardware, compute in FP32, at 3 to 6 times the GPU time. And
report a margin around any cut-off: here 150 to 334 of 3,893 variants sat within one median GPU difference of it,
and their calls can change with hardware alone. The forward test is necessary but not sufficient, because it passed
on every GPU.

## Limits

- The cut-offs were fitted and applied on the same assay data. They are research cut-offs for counting flips, not
  clinical thresholds. Calibrating a score for use as clinical evidence is a separate step
  [@pejaver2022calibration], and it was not done here.
- One gene, one model size (7B), one window length (8,192 bp) and one scoring rule. Flip rates depend on how many
  variants sit near the cut-off, so they will differ for other genes, cut-offs and score distributions.
- Identity within a GPU model was shown for the H100 only, on two machines at two providers. The L40S reproduced
  itself across hosts and drivers, but no second L40S provider was tested.
- Only three GPU products were tested. The other Evo 2 checkpoints (1B, 20B, 40B) need FP8 on Hopper GPUs
  and were not tested.
- The follow-up used 487 variants and one machine per GPU. Its subset cut-offs rest on 104 LOF variants, and it was
  exploratory.
- The cause was inferred from the pattern of results, not traced to kernels.
- Exact intervals treat variants as independent and the cut-off as fixed. SNVs at one position share a reference
  window, and flips cluster slightly (189 flips at 165 positions at the 95% cut-off), so the intervals may be a
  little narrow. The cut-offs were also fitted on the same scores.
- The AUROC difference interval was computed after the methods review, not pre-specified, and grouped bootstrap
  intervals may be too narrow because neighbouring positions share most of their window.

# Conclusion

Evo 2 7B gives the same scores on the same GPU model and different scores on a different one. The differences are
about $10^{-4}$. They moved AUROC by only 0.0021 (0.0002 to 0.0042) but changed 1.8% to 4.9% of BRCA1 calls at a
fixed cut-off. An L40 and an L40S disagree as much as an H100 and an L40S. FP32 reduces the difference about 250-fold,
at 3 to 6 times the cost. Validate Evo 2 on the GPU product that will run it, or compute in FP32.

# Data and code availability

Scores, receipts, logs, protocol and analysis code are in `rewire-bio/evo2-cloud-operations` at commit
`f8ae11c9a1c4` (<https://github.com/rewire-bio/evo2-cloud-operations>):

- `evidence/runs/recorded-b75c67aa/` and `evidence/runs/reproduction-30e7ae49/`: scores and receipts;
- `evidence/followup-20261011/`: the follow-up;
- `evidence/classification-flips/results.json`: the flip counts;
- `followup/classification_flips.py`: the script that makes them;
- `evidence/checks/verification.json`: the AUROC difference and top-ranked overlap.

The charts and this paper are in `rewire-bio/amp-2026-studies` (<https://github.com/rewire-bio/amp-2026-studies>),
folder `gpu-drift/`, in the first release tag after `amp-2026-v5` that contains `gpu-drift/paper/`. Study code is under the MIT licence.

Inputs keep their own licences. The Findlay et al. BRCA1 function scores are free for nonprofit use only; this is a
nonprofit use and the table is not redistributed [@findlay2018accurate]. Evo 2 code and weights are Apache-2.0
[@arcinstitute2026evo2]. Arc's test sequences and sample positions are Apache-2.0. The GRCh37 reference is public
domain.

