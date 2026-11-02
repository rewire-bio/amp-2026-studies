---
title: "Do genome foundation model scores help rank diagnostic variants in simulated rare-disease cases?"
subtitle: "The first run of a case benchmark: SNVs, indels, an inherited cancer panel and a 500-case replication"
author: '`Tim Richardson \\ Rewire`{=latex}'
date: 2 November 2026
abstract: |
  A clinical lab has to decide which of a patient's rare variants to look at first. We built simulated cases from
  public data to ask whether scores from genome foundation models help with that ranking. Each case is one 1000
  Genomes person's rare variants in the genes of a PanelApp panel, with one ClinVar pathogenic or likely pathogenic
  variant inserted. We ranked every variant with four established lookups (CADD, SpliceAI, REVEL, AlphaMissense),
  with four foundation-model tables (GPN-MSA, GPN-Star vertebrate and mammal, AlphaGenome AVI) and with Evo 2 7B,
  which we ran ourselves on one H100. We ran four studies: 200 SNV cases and 120 indel cases on the intellectual
  disability panel, 140 cases on an inherited cancer panel, and 500 new SNV cases from 50 new people. The four
  lookups combined in a logistic regression put the diagnostic variant first in 59.5% (95% interval 48.0 to 70.0)
  of the 200 SNV cases and in the top 5 in 92.0% (88.5 to 95.5). No foundation-model score alone did better.
  GPN-MSA's top-1 gain on the 200 cases (+6.0 points, +2.0 to +10.5) did not replicate on the 500 (+1.2, -3.0 to
  +5.8). AlphaGenome AVI's locked gain on the 500 was similar in size (+6.2, 0.0 to +13.0), but its interval reaches
  zero. It excluded zero only after a gnomAD AF 0.1% filter (+7.4, +1.2 to +14.2), not after an AF 1% filter (+5.2,
  -0.6 to +11.4) or a per-ancestry (grpmax) filter (+0.2, -3.6 to +3.8). Evo 2 added +15.8 points (+5.0 to +26.8)
  at top 1 on indels; the second run, with shorter lists and a higher baseline, does not reproduce this (+1.7, -5.9
  to +9.9). Review then found that the background people are in gnomAD, which flatters any frequency signal. A
  second run that fixes this finds no gain that matters from AVI (+1.8, -0.4 to +4.0). This paper reports the first
  run as it was and states which of its results the second run replaces.
---

# Introduction

A lab sequencing a child with intellectual disability finds hundreds of rare variants in the genes of the relevant
gene panel. It has to decide which to look at first. Labs filter on population frequency, use the predicted
consequence, and sort by a few established scores: CADD [@schubach2024cadd], SpliceAI [@jaganathan2019predicting],
REVEL [@ioannidis2016revel] and AlphaMissense [@cheng2023accurate]. Ranking tools such as Exomiser combine these
[@smedley2015next]. Genome foundation models now publish genome-wide score tables that claim to do better: GPN-MSA
[@benegas2025dna], GPN-Star [@ye2026predicting] and the AlphaGenome variant impact score, AVI
[@cheng2026alphagenome]. Evo 2 has no table, but a lab can run it [@brixi2026genome].

Most published comparisons of these scores measure how well they separate pathogenic from benign variants across
many genes. A lab's question is different. It has one patient and one list, and it wants the diagnostic variant
near the top of that list. This study measures that directly. It asks three questions:

1. How often does each score, alone, rank the diagnostic variant first or in the top 5?
2. How often do the established lookups, combined, do so?
3. Does adding one foundation-model score to the combined lookups raise those rates?

The cases are simulated from public data, so the right answer is known and no patient data are needed. This is the
first run of the benchmark, released at tag `amp-2026-v4` [@rewire2026studies]. After it was released, a review
found problems in the design (section 4.1). A second run fixes them [@rewire2026rdv2]. This paper reports the
first run's results unchanged and says which of them the second run replaces.

# Methods

The scripts in the repository directory `rare-disease/` are the full record. Each has a docstring with its
inputs, outputs and options. This section follows them.

## Design

Four studies share one design (Figure 1). A case is one 1000 Genomes person's rare variants in a panel's genes (the
background) plus one ClinVar pathogenic or likely pathogenic (P/LP) variant (the diagnostic variant). The diagnostic
variant is inserted into the list; the person does not carry it. Each score ranks every variant in the case, and the
outcome is where the diagnostic variant lands (Figure 2).

![The four studies. The SNV, indel and cancer studies use the same 50 people. Their combined models are fitted with
cross-validation. For the replication, models are fitted on all 200 SNV cases, locked, and applied unchanged to 500
cases from 50 new people. The follow-up checks were added after the first results.](figures/diagram_design.pdf){width=100%}

![How one case is built and scored.](figures/diagram_case_build.pdf){width=100%}

## Inputs

All coordinates are GRCh38. Versions are those used for the results.

- **Panels.** PanelApp panel 285, Intellectual disability, version 11.45, green genes only
  [@martin2019panelapp; @genomicsengland2026panelapp285]. It has 1,513 green genes, 1,512 of them with a MANE Select
  transcript. The cancer study uses panel 245, Adult solid tumours cancer susceptibility, version 2.40, with 82
  green genes [@genomicsengland2026panelapp245]. These versions can be fetched from the PanelApp API:
  `panelapp.genomicsengland.co.uk/api/v1/panels/285/?version=11.45` and `.../panels/245/?version=2.40`.
- **Regions.** The exons of each green gene's MANE v1.5 Select transcript [@morales2022joint], plus 20 bp on each
  side, merged (`regions.py`).
- **Backgrounds.** The 1000 Genomes high-coverage phased panel of 3,202 samples (release 20220422)
  [@byrskabishop2022high], fetched for the regions with bcftools [@danecek2021twelve] (`fetch_kgp.sh`). Fifty
  samples were drawn at random with seed 20261016 (`prepare_data.sh`). By superpopulation they are AFR 12, EUR 11,
  EAS 11, SAS 10 and AMR 6; 16 of them are among the 698 related samples. The SNV, indel and cancer studies all use
  these 50.
- **Diagnostic variants.** ClinVar, fileDate 2026-10-11 [@ncbi2026clinvar; @landrum2018clinvar]. A candidate is
  classified pathogenic, likely pathogenic or both, with at least two review stars (criteria provided by multiple
  submitters with no conflicts, expert panel, or practice guideline), and is at most 50 bp (`clinvar_plp.py`). Its
  consequence is the most severe term in ClinVar's MC field across transcripts, in VEP's order. `prepare_data.sh`
  downloads the current ClinVar and PanelApp files, not these versions. The ClinVar release used is the weekly
  archive file `vcf_GRCh38/archive_2.0/2026/clinvar_20261011.vcf.gz` on the NCBI FTP site. The committed
  `cases.tsv` files fix every case, so reproducing the results does not need these downloads.
- **Cancer inputs.** These were made by hand; `prepare_data.sh` does not cover them. From `data/cancer/`: download
  PanelApp panel 245 version 2.40 as `panelapp_245.json`; run `regions.py 245 cancer`; run the two `bcftools query`
  commands and `clinvar_plp.py` from `prepare_data.sh` step 3 against `cancer/regions_nochr.tsv`; concatenate the
  SNV and indel outputs into `cancer/clinvar_all_plp2.tsv`; and run `fetch_kgp.sh` with the same `samples50.txt`
  into `cancer/kgp/`.

## Building the cases

`build_cases.py` builds each study. A background is every variant of one person in the regions with a 1000 Genomes
allele frequency below 1%. Background variants that are themselves ClinVar P/LP with two or more stars are removed,
so each case has one right answer. Diagnostic variants are drawn at random to a quota per consequence group (Table
1). Draws use Python's `random` module with seed 20261016 (SNV), 20261017 (indel), 20261018 (cancer) and 20261023
(replication); after shuffling, case i is assigned to the (i mod 50)th person in sorted sample-ID order. Inheritance and phenotype are ignored: a heterozygous diagnostic variant in a recessive gene counts as a
diagnosis.

Table: The four studies. Groups are the diagnostic variant's consequence. People are 1000 Genomes samples; median
list is the median number of variants per case, counting the diagnostic one.

| Study | Cases | People | Diagnostic groups | Median list |
|---|---|---|---|---|
| SNV, panel 285 | 200 (4 per person) | 50 | 60 missense, 40 canonical splice, 30 nonsense, 70 synonymous, UTR or intronic | 234.5 |
| Indel, panel 285 | 120 | same 50 | 40 frameshift, 30 in-frame, 30 canonical splice, 20 UTR or intronic | 37 |
| Inherited cancer, panel 245 | 140 | same 50 | 40 missense, 25 each of nonsense, frameshift, canonical splice, non-coding | 17 |
| Replication, panel 285 | 500 (10 per person) | 50 new | 150 missense, 100 canonical splice, 75 nonsense, 175 synonymous, UTR or intronic | 208.5 |

**Indels.** Indels are ranked among the background indels in the same person only. The indel and cancer studies
oversample the rarer consequence groups so that each group has enough cases. Of the 140 cancer diagnostics, 30 are
indels. In the cancer study, the non-coding group also takes synonymous SNVs.

**Replication.** The replication was built on 23 October 2026, after the first results, to test them on new data.
Its 500 diagnostic SNVs are drawn with seed 20261023 from the same ClinVar pool, none at a position used as a
diagnostic in the 200, with the same consequence mix scaled to 500. Its 50 people are drawn with seed 20261023 from
the 2,504 unrelated samples, leaving out the first 50 and their families (2,425 eligible). By superpopulation they
are EAS 13, EUR 12, AMR 9, AFR 8 and SAS 8. Their genotypes were fetched in the same way, except chrX and chr15,
which were refetched in 8 parallel chunks after the first fetch failed (`REPLICATION.md`).

## Scores

Table 2 lists the scores. Each is looked up for every distinct case variant by `lookup_scores.py`, except Evo 2.
For every score, higher means more damaging; GPN-MSA, GPN-Star and Evo 2 have their sign flipped.

Table: The scores. "Lookups" are the four established scores labs already use.

| Score | Kind | Source | Covers |
|---|---|---|---|
| CADD v1.7 PHRED | Lookup | inclAnno files by remote tabix [@schubach2024cadd; @li2011tabix] | All SNVs; indels seen in gnomAD |
| SpliceAI | Lookup | Largest of the four delta scores, from CADD's annotation; SpliceAI 1.3.1 run by us on indels [@jaganathan2019predicting] | SNVs, indels |
| REVEL v1.3 | Lookup | `revel-v1.3_all_chromosomes.zip` from the REVEL site, hg38 column, highest over transcripts [@ioannidis2016revel] | Missense |
| AlphaMissense | Lookup | hg38 table, highest over transcripts [@cheng2023accurate] | Missense |
| GPN-MSA | Foundation-model table | Song Lab all-SNV table [@benegas2025dna; @songlab2025gpnmsa] | All SNVs |
| GPN-Star V and M | Foundation-model table | Song Lab tables, models v100-200m (vertebrate) and m447-200m (mammal), calibrated LLR [@ye2026predicting; @songlab2026gpnstar] | All SNVs |
| AlphaGenome AVI | Foundation-model table | Static all-SNV table, raw score [@cheng2026alphagenome; @googledeepmind2026downloads] | All SNVs |
| Evo 2 7B | Foundation model, run by us | One H100, section 2.4 [@brixi2026genome; @rewire2026evo2ops] | Every case variant except the 500 replication cases |

SpliceAI 1.3.1 was run with its bundled grch38 annotation, the default 50 bp distance and no masking (TensorFlow
2.15), on UCSC hg38 (`run_spliceai.sh`). The settings behind the SNV scores in CADD's annotation are not stated in
the CADD files, so the SNV and indel SpliceAI scores may not be on the same footing.

AVI was added to the SNV, replication and cancer studies on 24 October 2026, after the replication's first
results. It covers no indels, so the indel study has no AVI result, and AVI ranks every cancer indel diagnostic
last.

**Evo 2.** For each variant, two 8,192 bp windows of GRCh38, uppercased, with the variant at index 4,096: a reference
window and a variant window (`evo2/make_windows.py`). For an indel, the variant window keeps the reference length,
with the right flank refilled from the genome. Evo 2 7B gives the mean log-likelihood of each window, and the score
is minus the difference, variant minus reference, as in Arc Institute's BRCA1 notebook [@arcinstitute2026evo2].
The checkpoint is `evo2_7b`, and the windows come from the UCSC hg38 FASTA. Numerics were the defaults: bfloat16
weights, FP8 input projections, flash attention, batch size 1. Scoring ran on one rented NVIDIA H100 80GB HBM3 at
USD 3.99 per hour, in the container `ghcr.io/rewire-bio/evo2-cloud-operations` (digest in section 6), with scoring
code from `evo2-cloud-operations` commit `f8ae11c` (full hash in section 6) [@rewire2026evo2ops]. Three jobs scored
13,345, 198 and 2,814 variants (16,357 in all) in about 3 hours 31 minutes of pod time on two pods (the 198 were
scored on the first pod), for about USD 14, at 2.9 windows per second and 16.8 GiB peak GPU memory
(`evo2/README.md`, `evo2/launch*.log`, `evo2/run*/receipt.json`). The 13,345 were scored for an earlier build of the
200 SNV cases; the 198 are the variants the final build added; the 2,814 are the indel and cancer case variants
not scored before. The replication cases were not scored with Evo 2.

## Ranking and outcomes

For each case and score, the competitors are the background variants that score equal to or higher than the
diagnostic variant. Ties count against the score. If the diagnostic variant has no score, every background variant
is a competitor, as a failed lookup would leave it at the bottom of a sorted list. A missing background score is
dropped. The outcomes are top 1, top 5 and top 10: the share of cases with fewer than 1, 5 or 10 competitors.

## Combined models

"Lookups combined" is a logistic regression (scikit-learn defaults: L2 penalty, C = 1, lbfgs solver; standardised
inputs; balanced class weights; up to 2,000 iterations) [@scikitlearn2026linear] over every variant in every case, with the diagnostic label as outcome. Its
inputs are CADD, SpliceAI, REVEL and AlphaMissense (a missing score set to 0, with a missing flag for SNVs), the
variant's consequence class, and every score-by-class product, so that a score can count for missense variants and
not for intronic ones. The class comes from CADD's consequence for SNVs and from `annotate_indels.py` (MANE Select)
for indels. For indels, CADD and precomputed SpliceAI depend on gnomAD presence, so the models treat them as missing
for every indel and use the SpliceAI scores we ran ourselves. "Lookups + X" adds one foundation-model score X to the
same inputs.

The models are fitted with 5-fold cross-validation grouped by person (GroupKFold, without shuffling) [@scikitlearn2026groupkfold], so
no case is scored by a model that saw any case from the same person. For the replication, each model is also
fitted once on all 200 SNV cases and applied unchanged to the 500 ("locked"). The locked models are the main
replication result; the cross-validated models refit on the 500 alone.

## Statistics

Intervals are 95% percentile bootstrap intervals over people: 1,000 resamples of the 50 people, each with all their
cases (`case_benchmark.py`, seed 20261016). A gain is the paired difference in top-k rate between "Lookups + X" and
"Lookups combined" on the same cases. Each rate and each gain uses its own 1,000 resamples; within a gain, both
models are evaluated on the same resampled people. The models are not refitted within the bootstrap, so the
intervals reflect sampling of the evaluated people only, not uncertainty in the fitted models; for locked models,
uncertainty from the 200 training cases is not included. We call a result statistically
significant at the 5% level when its 95% interval excludes zero. An interval that ends at exactly 0.0 does not
exclude zero. No interval is adjusted for multiple comparisons.

No protocol was written in advance. The replication was designed after the first results, to test one stated
result: GPN-MSA's top-1 gain on the 200 SNV cases. That is the one comparison planned before its data were seen.
AVI was added after the replication had been built and is exploratory. So are every per-type result, the indel and
cancer results, and the follow-up checks below. With nine scores, six combined models (and their locked versions
on the 500), three cut-offs and up to six groups per study, some intervals will exclude zero by chance.

## Follow-up checks

Two checks were added after the first results (`REPLICATION.md`).

**AVI training chromosomes.** AVI's training labels are gnomAD frequency classes on chromosomes 1, 4, 7, 8, 10, 13
and 15 [@cheng2026alphagenome], and the 1000 Genomes samples are in gnomAD [@koenig2024harmonized]. `avi_heldout.py`
repeats the single-score ranking within the training chromosomes and within the others: a case is assigned by its
diagnostic variant's chromosome and ranked only among background variants in the same set. It reports point
estimates without intervals.

**gnomAD frequency filter.** A diagnostic pipeline removes common variants before ranking [@whiffin2017using].
`gnomad_af.py` looks up each variant's gnomAD genome allele frequency (AF) and grpmax AF from MyVariant.info's
`gnomad_genome` field [@xin2016myvariant], which holds gnomAD v3.1.2 genomes [@gnomad2021v312]. Grpmax is the
highest AF over the afr, amr, eas, nfe and sas groups, the set gnomAD uses for genome grpmax. A variant that is not
found gets AF 0 and is kept. `case_benchmark.py --gnomad-max-af` drops every variant at or above the threshold,
drops cases whose diagnostic variant is removed, and refits every model, including the locked ones on the filtered
200. Three filters were run: AF below 1%, AF below 0.1% (the threshold the AlphaGenome Atlas preprint uses for its
own rare-disease test [@cheng2026alphagenome]) and grpmax AF below 0.1%. A first version of the grpmax filter also
included the Middle Eastern group; it was corrected on 1 November 2026, and the results here are the corrected ones.

## Software and reproducibility

`prepare_data.sh`, `build_cases.py`, `lookup_scores.py`, `annotate_indels.py`, `run_spliceai.sh`, `gnomad_af.py`,
`case_benchmark.py` and `avi_heldout.py` run in that order (Figure 3), with the cancer inputs made by hand as in
section 2.2. The analysis runs on a CPU. Package versions
were not pinned: the analysis was run with `uv run --with pandas --with numpy --with scikit-learn`. On 24 October
2026 a rerun reproduced `results_indel.json`, `results_cancer.json` and `results_replication.json` byte for byte. In
`results.json` every row reproduced except "Lookups + Evo 2" at top 1, which gave 63.0% against the 62.5% committed
before; scikit-learn 1.7.2 and 1.9.1 gave the same 63.0%, and the cause is not known. Adding AVI that day changed
the order of the bootstrap draws, because all intervals share one random stream, so most intervals moved, typically
by 0.5 points or less. The unfiltered results are from that rerun (`results/compare_amp-2026-v1.txt` lists every
change); the frequency-filtered results are from later runs of the same code (grpmax rerun 1 November 2026).

![The scripts, in the order they run.](figures/diagram_pipeline.pdf){width=100%}

# Results

Every number below is from `rare-disease/results/` at tag `amp-2026-v4`, except where another file is named. Rates are percentages with 95% intervals in
brackets; gains are in percentage points.

## Single scores

On the 200 SNV cases, AVI was the best single score: 49.0% (41.0 to 57.5) at top 1 and 87.0% (82.0 to 91.5) at top
5 (Table 3). The best established lookup alone was SpliceAI at top 1, 39.5% (31.5 to 47.5), and at top 5, 54.5%
(47.5 to 61.5). REVEL and AlphaMissense score almost only missense variants (65 and 58 of the 200 diagnostics), so
they leave nearly every other diagnostic at the bottom of the list. Evo 2 was among the weakest at top 1, 16.5%
(11.0 to 22.0). On the 500 replication cases AVI was again the best single score: 53.2% (46.2 to 60.8) at top 1
and 82.8% (79.6 to 86.0) at top 5.

```{=latex}
\begingroup\small
```

Table: Single scores and combined models, 200 SNV cases (`results.json`, `strata.all.scores`), with top 1 and top 5
on the 500 replication cases (`results_replication.json`; combined models locked). Evo 2 was not run on the 500.

| Score | 200: top 1 | 200: top 5 | 200: top 10 | 500: top 1 | 500: top 5 |
|---|---|---|---|---|---|
| CADD v1.7 | 28.0 (22.0 to 35.0) | 51.5 (45.0 to 58.0) | 63.5 (57.0 to 69.5) | 30.2 (24.8 to 36.0) | 55.8 (50.4 to 61.0) |
| SpliceAI | 39.5 (31.5 to 47.5) | 54.5 (47.5 to 61.5) | 57.0 (50.0 to 64.0) | 40.2 (35.2 to 45.2) | 53.6 (49.0 to 58.4) |
| REVEL | 17.5 (12.5 to 22.5) | 27.0 (21.5 to 33.0) | 30.5 (25.0 to 37.0) | 16.4 (13.2 to 20.0) | 28.2 (24.8 to 31.6) |
| AlphaMissense | 17.5 (12.0 to 23.5) | 26.5 (20.5 to 32.5) | 28.5 (22.5 to 34.0) | 20.2 (16.8 to 23.6) | 27.2 (23.4 to 31.0) |
| GPN-MSA | 19.0 (13.5 to 25.0) | 48.5 (42.0 to 55.5) | 64.5 (59.0 to 70.0) | 23.0 (19.2 to 27.2) | 50.6 (45.2 to 55.6) |
| GPN-Star-V | 30.0 (23.5 to 36.0) | 57.5 (53.0 to 62.0) | 67.0 (62.0 to 72.0) | 36.6 (30.8 to 42.0) | 63.8 (58.4 to 69.0) |
| GPN-Star-M | 28.5 (22.0 to 35.0) | 52.5 (46.5 to 58.5) | 65.0 (60.0 to 70.5) | 29.4 (24.2 to 34.2) | 57.8 (53.2 to 62.4) |
| AlphaGenome AVI | 49.0 (41.0 to 57.5) | 87.0 (82.0 to 91.5) | 93.5 (90.0 to 96.5) | 53.2 (46.2 to 60.8) | 82.8 (79.6 to 86.0) |
| Evo 2 7B | 16.5 (11.0 to 22.0) | 43.0 (36.0 to 50.0) | 54.5 (48.0 to 61.5) | | |
| Lookups combined | 59.5 (48.0 to 70.0) | 92.0 (88.5 to 95.5) | 95.5 (92.5 to 98.0) | 61.0 (53.6 to 67.2) | 91.0 (88.6 to 93.2) |
| Lookups + AVI | 67.0 (57.5 to 76.0) | 93.5 (90.0 to 96.5) | 96.0 (93.0 to 98.5) | 67.2 (60.6 to 73.8) | 93.8 (92.0 to 95.4) |

```{=latex}
\endgroup
```


## The lookups combined

The four lookups combined with consequence class did better than any single score. On the 200 SNV cases they put
the diagnostic variant first in 59.5% (48.0 to 70.0) of cases, in the top 5 in 92.0% (88.5 to 95.5) and in the top
10 in 95.5% (92.5 to 98.0). With the model locked on the 500 they reached 61.0% (53.6 to 67.2) at top 1 and 91.0%
(88.6 to 93.2) at top 5 (Table 3).

Within each consequence group the combined model never beat the one lookup suited to the group at top 1 (Figure 4). On the 200 cases it
lost to the one lookup suited to the group in three of four groups: canonical splice 82.5% (68.4 to 94.4) against
SpliceAI 90.0% (80.5 to 97.6); nonsense 60.0% (40.6 to 77.4) against CADD 90.0% (76.6 to 100.0); missense 43.3%
(27.9 to 60.0) against REVEL 58.3% (44.3 to 71.2) and AlphaMissense 58.3% (44.6 to 71.7). On synonymous, UTR or
intronic diagnostics it tied with SpliceAI at 60.0% (44.8 to 72.3, and 47.1 to 72.2). These per-group differences were not tested with paired intervals. The combined model ranks all
variant types on one scale, which costs it within each type.

```{=latex}
\begin{figure}[p]
\centering
\includegraphics[angle=90,height=0.9\textheight,width=\textwidth,keepaspectratio]{figures/rd_by_type.pdf}
\caption{Top 1 by consequence group, 200 SNV cases. The hatched bar is the one lookup suited to the group, alone. The
line under each title is the gain from adding AVI to the lookups, with its 95\% interval.}
\end{figure}
```

## Adding a foundation-model score to the lookups

On the 200 SNV cases, adding GPN-MSA raised top 1 by +6.0 points (+2.0 to +10.5), which is significant at the 5%
level. Adding AVI raised it by +7.5 (-0.5 to +15.5), GPN-Star-V by +4.0 (-1.5 to +11.0), GPN-Star-M by +1.5 (-2.0
to +5.5) and Evo 2 by +3.5 (-1.5 to +9.0); none of these intervals excludes zero. At top 5 and top 10 no added
score changed the rate by more than 2.5 points (Figure 5, Table 4). The only top-5 change whose interval excluded
zero was a loss: adding GPN-Star-V lowered top 5 by 2.5 points (-5.5 to -0.5).

![Top 1, 5 and 10 for the lookups combined, alone and with AVI, GPN-MSA or Evo 2 added, 200 SNV cases. Bars show
95% intervals.](figures/rd_combined.pdf){width=100%}

```{=latex}
\begingroup\small
```

Table: Gain from adding each score to the lookups combined, points (`added_to_lookups` and
`added_to_locked_lookups`). The 200 are cross-validated; the 500 are shown with the models locked and with them
refitted by cross-validation on the 500.

| Added score | 200, top 1 | 500 locked, top 1 | 500 refitted, top 1 | 500 locked, top 5 | 500 locked, top 10 |
|---|---|---|---|---|---|
| GPN-MSA | +6.0 (+2.0 to +10.5) | +1.2 (-3.0 to +5.8) | +2.4 (-1.2 to +6.0) | +0.6 (-1.0 to +2.6) | +0.2 (-1.0 to +1.6) |
| GPN-Star-V | +4.0 (-1.5 to +11.0) | +3.0 (-1.2 to +8.0) | +5.6 (+1.8 to +9.8) | +0.6 (-1.0 to +2.2) | +0.4 (-0.8 to +1.8) |
| GPN-Star-M | +1.5 (-2.0 to +5.5) | +2.6 (-2.0 to +7.2) | +6.8 (+2.2 to +11.4) | +0.8 (-0.6 to +2.2) | +2.0 (+1.0 to +3.2) |
| AlphaGenome AVI | +7.5 (-0.5 to +15.5) | +6.2 (0.0 to +13.0) | +10.4 (+4.2 to +17.4) | +2.8 (+0.8 to +5.0) | +2.0 (+0.4 to +3.6) |
| Evo 2 | +3.5 (-1.5 to +9.0) | | | | |

```{=latex}
\endgroup
```


## Replication on 500 new cases

GPN-MSA's gain did not replicate. With the model locked from the 200, it was +1.2 points (-3.0 to +5.8) at top 1 on
the 500; refitted on the 500, +2.4 (-1.2 to +6.0). Neither interval excludes zero. The original and replication
intervals overlap, and no test of their difference was run (Figure 6). This was the planned test. The original +6.0 was the largest of
four gains measured on 200 cases, so a smaller value on new cases is what selection alone would predict.

AVI's gain was +6.2 (0.0 to +13.0) locked and +10.4 (+4.2 to +17.4) refitted. The locked interval reaches zero, so
the locked gain is not significant at the 5% level; the refitted one is. With the model locked, AVI's gains at top
5, +2.8 (+0.8 to +5.0), and top 10, +2.0 (+0.4 to +3.6), excluded zero, as did GPN-Star-M's at top 10, +2.0 (+1.0
to +3.2). Refitted, AVI's gain at top 5, +2.2 (+0.8 to +3.6), and GPN-Star-V's at top 10, +1.0 (+0.2 to +2.0),
also excluded zero. All are under 3 points. GPN-Star-V and GPN-Star-M had significant refitted gains at top 1,
+5.6 (+1.8 to +9.8) and +6.8 (+2.2 to +11.4), but not locked ones, +3.0 (-1.2 to +8.0) and +2.6 (-2.0 to +7.2).
All AVI and GPN-Star results here are exploratory.

![Gain at top 1 from adding each score to the lookups: the first 200 cases (hollow dots, cross-validated) and the 500
new cases with the model locked (filled dots). Lines are 95% intervals.](figures/rd_gains.pdf){width=100%}

On missense diagnostics the combined model again trailed AlphaMissense alone: on the 500, 35.3% (24.6 to 46.7)
refitted and 48.0% (37.7 to 57.6) locked, against 67.3% (58.8 to 75.0).

## AVI and its training chromosomes

AVI's lead did not depend on its training chromosomes (Table 5). On the other chromosomes, AVI put the diagnostic
variant in the top 5 in 91.1% of 157 cases from the 200 and 88.1% of 387 cases from the 500, ahead of the next
single score in both. Within the training chromosomes most scores did better, probably because the lists are
shorter (median 75 and 63 variants, against 163 and 149); SpliceAI on the 500 and GPN-MSA on the 200 did not. `avi_heldout.json` has no intervals, so these are point estimates.

Table: Top 5 for single scores, by chromosome set (`avi_heldout.json`). Lists are ranked within each set only, so
compare scores within a column, not across columns.

| Score | 200: training (43 cases) | 200: other (157) | 500: training (113) | 500: other (387) |
|---|---|---|---|---|
| CADD v1.7 | 65.1 | 57.3 | 77.0 | 58.9 |
| SpliceAI | 65.1 | 56.1 | 50.4 | 57.9 |
| GPN-MSA | 60.5 | 60.5 | 67.3 | 58.9 |
| GPN-Star-V | 72.1 | 61.1 | 76.1 | 68.7 |
| AlphaGenome AVI | 97.7 | 91.1 | 93.8 | 88.1 |

## gnomAD frequency filters

AF below 0.1% removed about 60% of each background (median list 234.5 to 123.5 on the 200, and 208.5 to 120.5 on
the 500) and no diagnostic variant. AVI's gain stayed similar in size, and its locked interval on the 500 moved
above zero: +7.4
(+1.2 to +14.2). Under the AF 1% filter it was +5.2 (-0.6 to +11.4). GPN-MSA's gain on the 200 fell from +6.0 to
+2.0 (-3.0 to +7.0), and to -1.0 (-8.0 to +5.0) under the AF 1% filter. This suggests, but does not show, that part
of it came from common background variants (Table 6).

The grpmax filter, below 0.1% in each of the five main ancestry groups, is the one a published pipeline applies:
the Genomics England tiering and Exomiser compare frequencies per population [@genomicsengland2025tiering;
@smedley2015next]. It left a median list of 51 on the 200 and 34.5 on the 500, and removed the cases whose
diagnostic variant was itself above the threshold (one of the 200, two of the 500). The lookups combined then put
the diagnostic variant first in 74.4% (63.0 to 83.5) of the 200 and 78.3% (71.2 to 84.4) of the 500 locked. No
added score had a top-1 gain whose interval excluded zero on either set. AVI's locked gain on the 500 was +0.2 (-3.6
to +3.8) (Figure 7).

![AVI added to the lookups, gain at top 1, as the tests get stricter. Filled dots are the 500 new cases with the
model locked. Lines are 95% intervals.](figures/rd_stricter.pdf){width=100%}

```{=latex}
\begingroup\small
```

Table: Gain at top 1 from adding a score to the lookups combined under each frequency filter, points. "Lookups"
is the lookups combined alone at top 1. Locked on the 500. Files: `results[_replication]_gnomad<af>[_grpmax].json`.

| Cases, filter | Median list | Lookups | + AVI | + GPN-MSA | + Evo 2 |
|---|---|---|---|---|---|
| 200, none | 234.5 | 59.5 (48.0 to 70.0) | +7.5 (-0.5 to +15.5) | +6.0 (+2.0 to +10.5) | +3.5 (-1.5 to +9.0) |
| 200, AF below 1% | 217 | 57.0 (46.0 to 67.0) | +9.5 (+1.0 to +18.0) | -1.0 (-8.0 to +5.0) | +1.5 (-4.0 to +7.5) |
| 200, AF below 0.1% | 123.5 | 64.5 (52.5 to 75.0) | +8.0 (+0.5 to +16.0) | +2.0 (-3.0 to +7.0) | +1.5 (-3.5 to +7.0) |
| 199, grpmax below 0.1% | 51 | 74.4 (63.0 to 83.5) | +4.0 (-3.0 to +12.1) | -1.0 (-7.0 to +4.6) | +1.0 (-4.5 to +6.5) |
| 500, none | 208.5 | 61.0 (53.6 to 67.2) | +6.2 (0.0 to +13.0) | +1.2 (-3.0 to +5.8) | |
| 500, AF below 1% | 200 | 61.4 (54.4 to 67.8) | +5.2 (-0.6 to +11.4) | +1.4 (-2.8 to +5.6) | |
| 500, AF below 0.1% | 120.5 | 66.2 (58.2 to 73.8) | +7.4 (+1.2 to +14.2) | +1.2 (-3.6 to +6.0) | |
| 498, grpmax below 0.1% | 34.5 | 78.3 (71.2 to 84.4) | +0.2 (-3.6 to +3.8) | -0.2 (-3.4 to +2.8) | |

```{=latex}
\endgroup
```


## Indels

On the 120 indel cases, ranked among a median of 37 indels, Evo 2 was the best single score: 57.5% (47.5 to 67.2) at
top 1 and 92.5% (86.8 to 96.8) at top 5. SpliceAI, run ourselves, reached 41.7% (32.5 to 51.2) and 62.5% (54.0 to
70.8); CADD scored only 36 of the 120 diagnostic indels and reached 17.5% (10.8 to 25.4) at top 1. The table-based
foundation-model scores cover no indels. The lookups combined reached 61.7% (51.2 to 71.3) at top 1 and 95.0% (90.1
to 99.1) at top 5. Adding Evo 2 raised top 1 to 77.5% (67.2 to 85.5), a gain of +15.8 points (+5.0 to +26.8), which
is significant at the 5% level (Figure 8). The gain at top 5 was +0.8 (-2.7 to +5.0).

![Top 1, 5 and 10 for the lookups combined, alone and with Evo 2 added, 120 indel cases. Bars show 95%
intervals.](figures/rd_indel_combined.pdf){width=100%}

By group, the gain at top 1 was in-frame +30.0 (+10.3 to +46.9, 30 cases), UTR or intronic +25.0 (+5.0 to +47.6,
20 cases), frameshift +12.5 (-5.4 to +31.8, 40 cases) and canonical splice 0.0 (-10.3 to +10.3, 30 cases). After the
gnomAD AF 0.1% filter (median list 15) the lookups reached 71.7% (61.0 to 80.7) and the gain halved to +8.3 (+0.9 to
+16.3), now mostly in frameshift +17.5 (+2.4 to +33.3), with in-frame +10.0 (-3.3 to +24.1) and the other two groups
at 0.0. With 20 to 40 cases per group, which group the gain comes from is not settled. All indel results are
exploratory.

## Inherited cancer panel

The cancer panel has 82 green genes, and its cases are short: a median of 17 variants. The lookups combined put the
diagnostic variant first in 85.0% (75.7 to 92.4) of the 140 cases, in the top 5 in 97.9% (95.1 to 100.0) and in the
top 10 in all of them. Evo 2 was the best single score at top 1, 72.9% (64.8 to 80.4), with AVI at 71.4% (62.5 to
79.6); AVI scores no indel and so ranks the 30 indel diagnostics last. CADD alone reached 56.4% (46.3 to 66.4).

Adding a score to the lookups raised top 1 by +2.9 to +5.0 points: GPN-MSA +2.9 (+0.7 to +5.8), GPN-Star-V +5.0
(+1.4 to +10.1), GPN-Star-M +3.6 (0.0 to +7.2), AVI +3.6 (0.0 to +7.3) and Evo 2 +4.3 (+0.7 to +9.3). The
GPN-MSA, GPN-Star-V and Evo 2 intervals exclude zero; the AVI and GPN-Star-M intervals end at zero. The gains were
largest on non-coding diagnostics, where the lookups reached 68.0% (48.0 to 87.5) and adding AVI gave +16.0 (+3.8 to
+30.8, 25 cases) (Figure 9). After the gnomAD AF 0.1% filter (median list 9) the lookups reached 87.9% (79.4 to
95.1) at top 1 and 99.3% (97.8 to 100.0) at top 5, and no added score gained more than 3 points: AVI +2.9 (0.0 to
+8.6), Evo 2 +2.9 (0.0 to +7.7). All cancer results are exploratory.

```{=latex}
\begin{figure}[p]
\centering
\includegraphics[angle=90,height=0.9\textheight,width=\textwidth,keepaspectratio]{figures/rd_cancer_by_type.pdf}
\caption{Top 1 by consequence group, 140 inherited cancer cases. AVI has no score for indels, so its frameshift bar
is 0.}
\end{figure}
```

**PP3 thresholds.** ClinGen's calibrated thresholds turn a missense score into PP3 evidence at Supporting, Moderate
or Strong [@pejaver2022calibration; @bergquist2025calibration; @richards2015standards]: AlphaMissense 0.792, 0.906
and 0.990, and REVEL 0.644, 0.773 and 0.932. Table 7 gives the share of missense variants reaching each level.
Background rows count a variant once in each case it appears in. The Supporting and Moderate thresholds were met
by most diagnostic variants and few background variants; Strong was met by 28 to 75% of diagnostic variants and at
most 0.4% of background ones. The diagnostic rates may be too high, because the diagnostic variants are
ClinVar P/LP classifications that may have used these tools as PP3 evidence, and REVEL was trained on pathogenic
missense variants [@grimm2015evaluation; @ioannidis2016revel]. The rates differ between panels, and no intervals
were computed.

Table: Share of missense variants at or above each PP3 threshold, percent (`pp3_missense` in each results file).
Missense is CADD's consequence class; a variant with no score counts as below every threshold (in the 200 cases, no
AlphaMissense score for 328 of 10,196 background rows and 2 of 60 diagnostic rows, and no REVEL score for 124
background rows).

| Panel, variants | n | AlphaMissense: Sup., Mod., Strong | REVEL: Sup., Mod., Strong |
|---|---|---|---|
| Cancer, diagnostic | 40 | 95.0, 90.0, 75.0 | 85.0, 70.0, 40.0 |
| Cancer, background | 483 | 3.3, 1.9, 0.0 | 6.8, 2.9, 0.4 |
| Intellectual disability (200), diagnostic | 60 | 78.3, 65.0, 33.3 | 88.3, 71.7, 28.3 |
| Intellectual disability (200), background | 10,196 | 3.5, 1.8, 0.3 | 5.3, 2.0, 0.2 |

# Discussion

Across four studies, the established lookups combined with consequence class did most of the work. No
foundation-model score on its own matched them at top 1 or top 5 on any SNV study. AVI was the best single score,
but it is not a single model output: it combines 18 features, including AlphaGenome's predicted effects,
AlphaMissense, conservation and consequence flags [@cheng2026alphagenome], so its lead over single scores is partly
a combination against single scores.

The one planned test, GPN-MSA's gain on new cases, failed: +6.0 on 200 cases became +1.2 on 500. That pattern, the
largest of several estimates shrinking on new data, applies to every exploratory gain in this paper. AVI's gain was
similar in size on the 500, but its locked interval excluded zero only under the AF 0.1% filter, and it was this
result that the talk built on. It vanished under the grpmax filter, where the lists are short and the lookups
already rank the diagnostic first in three cases in four. On indels, Evo 2 was the strongest single score and gave
the largest gain of any added score, +15.8 points, halving to +8.3 after the frequency filter; the second run does
not reproduce it (section 4.1). On the inherited cancer panel the lists are short and the
lookups leave little room for any added score.

## What the second run found

After release, review of the study found five problems, recorded as issues #1 and #2 in the repository. The
second run's proposal adds a sixth.

1. **The background people are in gnomAD** (issue #1). 1000 Genomes is part of gnomAD [@koenig2024harmonized], so
   all but 14 of the 13,145 distinct background SNVs in the 200 cases, and all 11,208 in the 500, had been seen
   there (`data/gnomad_af.tsv`), while 75% of the 200 diagnostic SNVs and 76% of the 500 had not. "Not in
   gnomAD" alone picked out the diagnostic variant in about three cases in four. A real patient's private variants
   are absent from gnomAD and compete with the diagnostic one. This flatters any score that has learned frequency.
   AVI was trained on gnomAD frequency classes. It also makes the gnomAD-filtered lists easier than a real
   patient's, so the lookups' 74 to 80% top 1 under grpmax may be too high.
2. **The realistic frequency filter was a side check** (issue #2). The study called the overall AF filter the one a
   pipeline would apply, but published pipelines filter per ancestry group, which is the grpmax filter. Under it,
   AVI's locked gain was +0.2 (-3.6 to +3.8).
3. **The combined model was a weak baseline within each variant type** (issue #2). It lost to the matching single
   tool in three of four groups (section 3.2), so part of any added score's gain could be making up for the model.
4. **Two consequence annotations were mixed** (issue #2). Groups came from ClinVar's MC field and model classes from
   CADD, so a variant could be in one group in the results and a different class in the model.
5. **Small groups and unpaired reporting** (issue #2). With 20 to 70 cases per group, most per-group differences
   are noise.
6. **Carrier states** (`rare-disease-v2/PROPOSAL.md`). The second run's proposal adds this problem: many diagnoses
   were not diagnoses. 94 of the 200 and 231 of the 500 diagnostic variants were single heterozygous variants in
   genes PanelApp lists as biallelic only, so they were carrier states.

The second run, in `rare-disease-v2/` at tag `amp-2026-v5` [@rewire2026rdv2], fixes these on the intellectual
disability panel only. It computes gnomAD frequencies without each person and their family, makes grpmax below 0.1%
the main filter, builds an inheritance filter into the cases, uses one consequence annotation, adds a rule-based
baseline, and fixes its analysis and reading rule before scoring its 500 new cases. Its primary result: adding AVI
to the lookups raises top 1 from 90.6% (87.8 to 93.4) to 92.4% (89.8 to 94.8), a gain of +1.8 points (-0.4 to
+4.0), which by its rule fixed in advance is no gain that matters. Its lists are short (median 19 variants) and the
lookups already rank the diagnostic first in nine cases in ten. GPN-MSA and GPN-Star add nothing. It agrees with
this run's grpmax result, not with its headline.

So this run's claim that AVI's gain holds on new cases should no longer be used. The second run replaces it. The
first run's results that the second run does not retest are the inherited cancer panel results and the cancer PP3
rates. They share problems 1, 3 and 4, and they are exploratory with small groups.

The second run does retest Evo 2 on indels, scored after its other results (exploratory). On its 120 indel cases,
with lists of a median 4 indels and a lookups baseline of 84.2% (75.4 to 92.4), adding Evo 2 gains +1.7 points at
top 1 (-5.9 to +9.9) (`rare-disease-v2/results/indel_main.json`, commit `6697db7`). The first run's +15.8 does not
hold under that design. That result was committed after tag `amp-2026-v5`, so it is cited by commit.

## Limits

- **gnomAD presence** (problem 1 above). Every frequency-related result here, including AVI's lead and the
  filtered results, is affected.
- **Inserted diagnoses.** The diagnostic variant is inserted, not observed. Inheritance, zygosity, phenotype and
  segregation are ignored, and many diagnoses were carrier states (problem 6).
- **Backgrounds are not patients.** 1000 Genomes participants were not selected for disease, so the backgrounds are
  assumed to hold no diagnosis. Only ClinVar P/LP variants with two or more stars were removed from them.
- **Circularity.** Many ClinVar P/LP classifications used computational evidence, and some scores were trained on
  ClinVar labels [@grimm2015evaluation]. This favours those scores on the diagnostic variants.
- **Few people.** Each study has 50 people, so intervals are wide. The replication uses 50 different people, with
  fewer African backgrounds and so shorter lists.
- **The combined model** is a plain logistic regression that loses to single tools within types (problem 3).
- **Multiplicity.** Apart from the GPN-MSA replication, every comparison is exploratory and unadjusted.
- **Software versions** were not pinned, and one row of `results.json` did not reproduce exactly (section 2.9).
- **Shared bootstrap stream.** Adding a score changed every interval slightly, because all intervals draw from one
  random stream.
- **Model uncertainty.** The bootstrap does not refit the models, so no interval includes the uncertainty from
  fitting them. This matters most for the locked replication gains, fitted on 200 cases from 50 people.

# Conclusion

In this first run, the established lookups combined did most of the work of ranking a simulated case's diagnostic
variant. GPN-MSA's early gain did not replicate. AVI's gain was similar in size on new cases, but its interval
excluded zero only under one of three frequency filters, and it vanished under the per-ancestry filter a pipeline
uses. Evo 2 helped most on indels, where few tables exist. The second run, with shorter lists and a higher
baseline, finds no indel gain that matters (+1.7, -5.9 to +9.9). The design flattered frequency, because the background people are in gnomAD. The second run, which fixes
that, finds no gain that matters from AVI, and its results replace this run's AVI claim.

# Data and code availability

Code, case files, per-case scores (except CADD and REVEL), Evo 2 scores and receipts, result files and logs are in
the repository `github.com/rewire-bio/amp-2026-studies` (MIT licence), directory `rare-disease/`, at tag
`amp-2026-v4` [@rewire2026studies]. The results come from `rare-disease/results/*.json` at that tag; Evo 2 run times
and costs come from `rare-disease/evo2/`, and gnomAD presence counts from `rare-disease/data/gnomad_af.tsv`.
`rare-disease/README.md` gives the order in which to run the scripts, and `rare-disease/data/README.md` describes
every data file. The second run is in `rare-disease-v2/` at tag `amp-2026-v5` [@rewire2026rdv2]; its Evo 2 indel
result is at commit `6697db7`.

The Evo 2 container is
`ghcr.io/rewire-bio/evo2-cloud-operations@sha256:`\
`22970720ff8bc05a00460a6066a70ed71951001aa52965141591a0ca6dcc968b`, and the scoring code is
`evo2-cloud-operations` commit `f8ae11c9a1c45505487c978b4cb7e064c4e59905` (`rare-disease/evo2/run*/instance.json`
and `ops_commit.txt`).

The inputs are public. CADD and REVEL are free for non-commercial use but not licensed for redistribution, so their
per-case scores are not committed and `lookup_scores.py` fetches them again. The AlphaMissense table is CC BY-NC-SA
4.0 and is downloaded, not committed. The AlphaGenome AVI table is downloaded under the terms on its downloads page
[@googledeepmind2026downloads]. PanelApp, ClinVar, MANE, the 1000 Genomes high-coverage panel, gnomAD v3.1.2 (via
MyVariant.info), and the GPN-MSA and GPN-Star tables are public. Evo 2 7B weights are from Arc Institute
[@arcinstitute2025evo2].
