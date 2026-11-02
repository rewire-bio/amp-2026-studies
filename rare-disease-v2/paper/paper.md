---
title: "Do genome foundation model scores help rank diagnostic variants once frequencies are realistic?"
subtitle: "A second, pre-specified run of a rare-disease case benchmark"
author: '`Tim Richardson \\ Rewire`{=latex}'
date: 2 November 2026
abstract: |
  The developers of genome foundation models report that their scores rank pathogenic variants better than
  established tools. A first run of this benchmark found a top-1 gain from AlphaGenome's variant impact score (AVI),
  but its background people were in gnomAD and about half its diagnoses were carrier states. This second run leaves
  each person and their family out of gnomAD, filters by grpmax below 0.1% and by mode of inheritance, and tests 500
  cases from 50 new 1000 Genomes people with models fitted on 200 earlier cases and locked. Adding AVI to the
  combined lookups (CADD, SpliceAI, REVEL, AlphaMissense) raised top 1 from 90.6% (95% interval 87.8 to 93.4) to
  92.4% (89.8 to 94.8), a gain of +1.8 points (-0.4 to +4.0): no gain that matters under the rule fixed in advance.
  Added to a rule-based ACMG points ranking, AVI gave +2.6 (+0.4 to +5.6; unadjusted, with no test pre-specified),
  mostly in non-coding diagnostics (+8.6, +3.2 to +15.1, descriptive). GPN-MSA, GPN-Star vertebrate and GPN-Star
  mammal added no gain that matters at top 1: every pre-specified 98.33% interval on the 500 spans zero and lies within -1.8 to
  +2.4. Evo 2, added after the other results were seen, gave +0.8 (98.75% interval -0.6 to +2.4) on the 500 and
  +1.7 (-9.0 to +12.8) on 120 indel cases; on the indel cases its top-5 gain, +4.2 (98.75%: +0.8 to +9.5), excludes
  zero. The first run reported +15.8 (95%: +5.0 to +26.8) on its indel cases, a different case set whose baseline
  lacked a CADD score for most indel diagnostics; the intervals overlap. With realistic filters, the established
  lookups rank the diagnostic first in about 9 cases in 10, which leaves added scores little room.
---

# Introduction

A clinical lab sequencing a child with intellectual disability finds tens of rare variants in the genes of the
relevant panel. It has to decide which to look at first. Labs rank them with a frequency filter, a consequence class
and a few established scores: CADD [@schubach2024cadd], SpliceAI [@jaganathan2019predicting], REVEL
[@ioannidis2016revel] and AlphaMissense [@cheng2023accurate], ranked by tools such as Exomiser [@smedley2015next].
Foundation models now publish genome-wide score tables, and their papers report better variant ranking: AlphaGenome's variant impact score (AVI)
[@cheng2026alphagenome], GPN-MSA [@benegas2025dna], GPN-Star [@ye2026predicting], Evo 2 [@brixi2026genome] and
EVEE [@pearce2026interpretable]. The question for a lab is whether adding one of them to what it already uses puts
the diagnostic variant first more often.

A first run of this benchmark [@rewire2026studies] built cases from a real 1000 Genomes person's rare variants in
the PanelApp intellectual disability panel, added one ClinVar pathogenic variant, and asked where each score ranked
it. AVI was the one added score whose top-1 gain stayed positive on 500 new cases: +6.2 points with models fitted on
the first 200 cases and locked, with an interval that reached zero (0.0 to +13.0). A review then found that the design favoured frequency, and any score that had
learned it:

- The background people are in gnomAD. 1000 Genomes is part of gnomAD, so every background variant had been seen
  there at least once, while about three quarters of the inserted ClinVar variants had not. "Not in gnomAD" alone
  picked out the diagnostic variant in most cases. AVI was trained on gnomAD frequency classes.
- About half the diagnoses were not diagnoses. 94 of the 200 and 231 of the 500 diagnostic variants were single
  heterozygous variants in genes that PanelApp lists as biallelic only: carrier states.
- The realistic filter, a maximum allele frequency per ancestry group (grpmax), was only a side check. Under it,
  AVI's locked gain on the 500 was +0.2 points (-3.6 to +3.8).
- The combined model was a weak baseline within each variant type. It lost to the matching single tool in three of
  four types on the 200 cases, so part of an added score's gain could be making up for that.
- Two consequence annotations were mixed: ClinVar's for the reported groups and CADD's for the model classes.

This second run fixes the last four of those problems and reduces the first. Leaving each person's family out of
gnomAD still leaves most diagnostics absent from gnomAD and almost every background variant present, so section 2.10
tests whether any gain depends on that. It asks four questions. Q1 (primary): with frequencies closer to what a
patient's pipeline would see, a grpmax filter below 0.1% and an inheritance filter, does adding AVI to the combined lookups raise top 1
on 500 cases from 50 people the first run never used? Q2: how often do the lookups alone rank the diagnostic variant
first, combined by a fitted model and by a rule-based ranking that follows ClinGen guidance? Q3: does AVI help the
rule-based ranking too? Q4 (secondary): do GPN-MSA, GPN-Star vertebrate or GPN-Star mammal add anything? The design,
the analysis and the rule for reading the primary result were written into a protocol before the replication cases
were scored; section 2.11 sets out which changes came after which results. The grpmax filter was made the main
filter knowing that AVI had shown no gain under it in the first run (+0.2 points). Evo 2 was added after the other
results had been seen and is reported as exploratory.

# Methods

The protocol (`PROTOCOL.md`) and the scripts in the repository are the full record. This section follows them. Every
change to the protocol after it was written is in its deviation log (section 2.11).

## Cases and people

A case is one 1000 Genomes person's rare variants in the panel's genes (the background) plus one ClinVar
pathogenic or likely pathogenic variant (the diagnostic variant), which is inserted and not carried by the person.
Each score ranks every variant in the case, and the outcome is where the diagnostic variant lands.

![Study design. Models are fitted on the 200 SNV cases from the first run's people, locked, and applied unchanged to
500 cases from 50 new people. The reading rule for the primary contrast was fixed before the 500 were scored.](figures/diagram_design.pdf){width=100%}

Figure 1 shows the design and Figure 2 how one case is built.

Inputs were fixed to the first run's. The panel is PanelApp panel 285, intellectual disability, version 11.45,
green genes only [@martin2019panelapp; @genomicsengland2026panelapp285]. The regions are the MANE v1.5 Select
[@morales2022joint] exons of those genes plus 20 bp on each side. Diagnostic candidates are ClinVar variants
(fileDate 2026-10-11) classified pathogenic or likely pathogenic with at least two review stars
[@ncbi2026clinvar]. Backgrounds come from the 1000 Genomes high-coverage phased panel (20220422)
[@byrskabishop2022high]. All coordinates are GRCh38.

**Background.** Each person's variants in the regions with 1000 Genomes allele frequency below 1%, then with grpmax
below 1% (section 2.2). Background variants that are ClinVar pathogenic or likely pathogenic at any review status
are removed, so each case has one right answer. The background keeps both SNVs and indels, because the inheritance
filter counts both (section 2.3); only the study's own variant type is ranked.

**Studies.** Three sets of cases (Table 1). The SNV and indel studies reuse the first run's 50 people. Indels are up
to 50 bp; variants with REF and ALT both longer than 1 bp are left out, as the 1000 Genomes backgrounds hold none
and SpliceAI cannot score them. The replication study uses 50 new people: the 2,504 unrelated 1000 Genomes samples,
minus both of the first run's sets of 50, minus anyone who shares a family ID with them or is their parent, child or
full sibling in the pedigree. The rest are sorted by ID, shuffled with `random.Random(20261110).shuffle`, and the
first 50 taken. Family IDs come from `20130606_g1k_3202_samples_ped_population.txt`, and the unrelated samples are
the 3,202 minus the 698 in `1000G_698_related_high_coverage.sequence.index`.

| Study | Cases | People | Per person | Diagnostic quotas | Seed |
|----------|-----|------------|--------|------------------------------|--------|
| SNV | 200 | first run's 50 | 4 | missense 60, canonical splice 40, nonsense 30, non-coding 70 | 20261102 |
| Indel | 120 | first run's 50 | 2 or 3 | frameshift 40, inframe 30, canonical splice 30, UTR or intronic 20 | 20261103 |
| Replication | 500 | 50 new people | 10 | missense 150, canonical splice 100, nonsense 75, non-coding 175 | 20261109 |

Table: The three studies. Non-coding SNV diagnostics are synonymous, UTR, intronic or splice region variants. The
replication cases are the primary test set.

**Diagnostic draw.** One Python `random.Random(seed)` with the study's seed is used for the groups in the order of
Table 1. Each group's eligible ClinVar rows (full gnomAD grpmax below 0.1%) are listed. For the replication study,
the list first drops every position used as a diagnostic in this run's 200 or in the first run's 200 or 500. The
list is sorted by chromosome (as text), position (as a number), REF and ALT, then permuted with
`rng.sample(list, len(list))`. The list is walked in that order, skipping a position already used in
the study and a gene that already has 10 diagnostics in that group, until the quota is met. The gene cap stops one
gene dominating a group (NF1 has 75 of the 515 eligible non-coding positions). Variants whose only ClinVar genes are
MSH3 or MUTYH are not eligible: they lie inside a panel gene's region (DHFR, TOE1) but cause a cancer predisposition outside the panel.
A group short of positions would have its quota reduced and the shortfall logged; none was. The drawn
diagnostics are then permuted with `rng.shuffle` and assigned round robin to the study's people, sorted by ID.

## Frequencies

A real patient is not in gnomAD, and usually neither are their relatives. Everyone else is. So for each case, the
gnomAD frequency of every background variant is recomputed without that person and their released relatives ("leave
the family out"). This is the main frequency.

![How one case is built. The diagnostic keeps its full gnomAD frequency, which is below 0.1% for every diagnostic, so
the frequency filter never removes it. It joins the list before the inheritance filter, so it counts toward the
two-variant rule.](figures/diagram_case_build.pdf){width=55%}

The source is gnomAD v3.1.2 genomes (76,156 genomes) [@chen2024genomic]. Full counts (AC and AN, overall and
for the groups afr, amr, eas, nfe and sas) come from MyVariant.info's `gnomad_genome` field. Here grpmax is the
highest allele frequency over those five groups, the groups gnomAD v4 uses for genome grpmax.

Family counts come from gnomAD's own genotype calls for the 1000 Genomes samples, in the HGDP and 1000 Genomes
subset release [@koenig2024harmonized] (`gnomad.genomes.v3.1.2.hgdp_tgp`, per-sample GT, GQ, DP and AD), read at each
person-variant position. Only samples that gnomAD released are counted, each in its inferred group. A genotype
counts only if it passes gnomAD's adj filter: GQ at least 20, DP at least 10 (haploid calls too), and for
heterozygotes an alternate allele balance (alternate AD over DP) of at least 0.2. The family is the person plus their
parents, children and full siblings in the 1000 Genomes pedigree. Relatedness not documented in the pedigree is not
handled.

For group *g*, the frequency is (AC~g~ minus AC~fam,g~) / (AN~g~ minus AN~fam,g~), with AC clamped at 0. AN~fam~ is
the adj-called alleles of the family members in that group: 2 per diploid call, 1 per haploid call on chrX outside
the pseudoautosomal regions. Where the subset file has no row, no family member carries the allele, so AC~fam~ is 0
and AN~fam~ is 2 (or 1) per released family member in the group. Family members in a group outside the five (one
replication person is in gnomAD's fin group) are subtracted from the overall counts only, so grpmax does not change
for them. A variant not in gnomAD, or with no allele left, has frequency 0 and is kept, as a pipeline keeps a variant
it cannot find. This includes 14 background variants of the first run's people that are in gnomAD's subset file
but have no MyVariant.info record (`results/gnomad_checks.md`, section 5); none of the replication people's variants
are affected. A gnomAD site that is not PASS is treated as present with its counts. Diagnostic variants keep their full gnomAD frequency, because the person does not carry them.

16 of the first run's 50 people are not in the gnomAD release (most of the 702 unreleased 1000 Genomes samples are
relatives removed by gnomAD's relatedness pruning). For them only released relatives are subtracted. All 50
replication people are in the release.

Checks before use (`results/gnomad_checks.md`): recounting, with the adj filter as above, the samples behind the
subset file's own AC and AN reproduces those counts on 32,445 of 32,445 rows (including 1,592 of 1,592 on chrX
outside the PARs). The subset file's full-gnomAD fields match MyVariant.info on 25,535 of 25,535 rows, overall and
in all five groups. The released 1000 Genomes counts never exceed gnomAD's group counts. 20 person variants (five on
chrX outside the PARs, one in a PAR) were worked by hand against the gnomAD sites VCF. No count needed clamping.

## Filters

**Frequency.** The main filter keeps a variant only if its grpmax is below 0.1%, in every gene. This follows the
dominant-mode threshold of Exomiser and of the Genomics England tiering, both applied per population
[@genomicsengland2025exomiser; @genomicsengland2025tiering]. Both are looser in other settings. Genomics England
tiering uses 1% in every group for recessive modes, and 0.2% (eas, sas) or 0.3% (asj) for some genome groups in
dominant modes; it also uses the gnomAD v2.1.1 exomes and the fin group [@genomicsengland2026popfreq]. Exomiser uses
2% for compound heterozygous recessive modes. The main filter is therefore stricter than Genomics England for the
277 diagnostics in biallelic genes (Table 9); the 1% check (section 2.10) matches its recessive threshold. The
filtering allele frequency and a disease-specific maximum credible allele frequency
[@whiffin2017using] were considered and not used. A looser background filter, grpmax below 1%, is a check (section
2.10). Cases are built once at 1% and every frequency is stored, so the main filter is applied at analysis time by
one function shared by every score. Every diagnostic has full grpmax below 0.1%, so no case loses its diagnostic
under either filter and the case counts are fixed before any score is looked at.

**Inheritance.** Each green gene's mode of inheritance comes from PanelApp, matched to the MANE gene of the region
(by Ensembl gene ID, then symbol, then alias). Each person's sex comes from the 1000 Genomes pedigree. Every case is
a diagnosis under its gene's PanelApp mode for intellectual disability (Table 2). The ClinVar condition and
inheritance asserted for each diagnostic variant were not checked against the panel, apart from excluding MSH3 and
MUTYH. Scores do not use zygosity; it decides only which background variants a
pipeline would keep.

| PanelApp mode | Diagnostic genotype | Background variant kept if |
|---|---|---|
| Monoallelic, both, or unknown (GNAS only) | heterozygous | always |
| X-linked, monoallelic in females | heterozygous (XX), hemizygous (XY) | always |
| X-linked, biallelic in females | homozygous (XX), hemizygous (XY) | XY: always; XX: as biallelic |
| Biallelic | homozygous | homozygous, or two or more variants in the gene |

Table: Diagnostic genotype and background inheritance filter by mode of inheritance. "Two or more variants" counts
the person's variants in the gene after the frequency filter, of both types and including the diagnostic; phase is
ignored, so this keeps possible compound heterozygotes. Variants outside every green gene's region are kept.

Homozygous diagnostics in biallelic genes fit a consanguineous family. Recessive diagnoses made of two different
variants are not modelled, because scores rank single variants.

## Consequence annotation

One annotation chooses the diagnostic variants, sets the fitted models' classes and defines the reported groups
(Table 3).
For SNVs it is CADD v1.7's most severe consequence over its rows for MANE v1.5 Select transcripts (CADD `FeatureID`,
version stripped); a variant with no MANE row would take the most severe over all rows, and none needed to. The
order, most severe first, is `STOP_GAINED`, `STOP_LOST`, `CANONICAL_SPLICE`, `NON_SYNONYMOUS`, `SPLICE_SITE`,
`SYNONYMOUS`, `5PRIME_UTR`, `3PRIME_UTR`, `NONCODING_CHANGE`, `INTRONIC`, `UPSTREAM`, `DOWNSTREAM`, `REGULATORY`,
`INTERGENIC`. For indels it is `annotate_indels.py`, which assigns the most severe class over overlapping MANE Select
transcripts, in this order: splice (overlaps the two intronic bases at an exon-intron boundary), frameshift (overlaps
coding sequence, length change not a multiple of 3), inframe (overlaps coding sequence, length change a multiple of
3), UTR (overlaps an exon outside the coding sequence), intronic (inside the transcript), other. It is run on every
indel, so diagnostic and background indels get the same information.

| Class | SNV (CADD consequence) | Indel (`annotate_indels.py`) |
|---|---|---|
| nonsense | `STOP_GAINED` | |
| frameshift | | frameshift |
| inframe | | inframe |
| missense | `NON_SYNONYMOUS` | |
| canonical splice | `CANONICAL_SPLICE` | splice |
| splice region | `SPLICE_SITE` | |
| synonymous | `SYNONYMOUS` | |
| UTR | `5PRIME_UTR`, `3PRIME_UTR` | utr |
| intronic | `INTRONIC` | intronic |
| other | `STOP_LOST` and anything else | other |

Table: Consequence classes. The SNV non-coding group is splice region, synonymous, UTR and intronic.

On the whole ClinVar SNV pool, CADD's class and ClinVar's consequence field disagree for 20 of about 22,000
canonical splice, nonsense and missense variants (`results/pool_counts.md`). The change of annotation from the
first run moves variants within the non-coding group (about a third become splice region) and moves ClinVar
initiator codon variants to missense.

For the fitted models only, synonymous, UTR and other are merged into one class, and so is any class with fewer than
10 diagnostics in the training cases. The SNV model classes are canonical splice, intronic, missense, nonsense,
splice region and other; the indel model classes are canonical splice, frameshift, inframe, intronic and other.

## Scores

Table 4 lists the scores.

| Score | Role | Variants scored | Source |
|---|---|---|---|
| CADD v1.7 PHRED [@schubach2024cadd] | lookup | SNVs | remote tabix, inclAnno table |
| SpliceAI [@jaganathan2019predicting] | lookup | SNVs (CADD's annotation); indels (run here) | max of the four deltas |
| REVEL v1.3 [@ioannidis2016revel] | lookup | missense SNVs | max over transcripts |
| AlphaMissense [@cheng2023accurate] | lookup | missense SNVs | hg38 table, max over transcripts |
| AlphaGenome AVI [@cheng2026alphagenome] | primary added score | SNVs | all-SNV table, raw score |
| GPN-MSA [@benegas2025dna] | secondary added score | SNVs | hg38 scores, minus LLR |
| GPN-Star V and M [@ye2026predicting] | secondary added scores | SNVs | vertebrate (v100) and mammal (m447) models, minus calibrated LLR |
| EVEE [@pearce2026interpretable; @pearce2026evee] | descriptive only | ClinVar variants | Zenodo 19701997, pathogenicity |
| Evo 2 7B [@brixi2026genome] | exploratory added score | SNVs and indels | scored here (below) |

Table: Scores. Higher means more damaging for every score as used.

A missing score ranks last, as a failed lookup would. CADD and precomputed SpliceAI are not used for indels in any
model or in the rule-based baseline: CADD's indel table has a row only for indels seen in gnomAD, which would give the
answer away. For the indel study's indels, SpliceAI 1.3.1 (`-A grch38`, default distance 50, unmasked; Python 3.11,
TensorFlow 2.15 and NumPy below 2, outside `requirements.txt`; the UCSC hg38 FASTA) scores every indel;
indels the first run had already scored with the same command and genome are reused (`run_spliceai.sh`). The score
is the largest of the four deltas over all annotated genes. EVEE's table has one row per ClinVar variant and is joined at EVEE position equal to VCF position minus
1 (its positions are 0-based). For every diagnostic with an EVEE row (200 in the 200, 492 in the 500, 69 in the
120), EVEE's variation ID equals the case's ClinVar ID. Because EVEE covers ClinVar variants only and its
released scores may include its training labels, it is analysed apart (section 2.7) and never compared with the main
analysis.

**Evo 2.** Evo 2 7B scores every case variant in all three studies, with the first run's method. For each variant, an
8,192 bp reference window is taken from hg38 with the variant at index 4,096, and a variant window is made from it.
Windows are uppercased, because the UCSC FASTA marks repeats in lower case and Evo 2's byte tokenizer treats the two
cases as different inputs. For an indel the variant window keeps the reference length: REF is replaced by ALT and the
right flank is refilled from the genome, so mean log-likelihoods compare like for like. Each REF is checked against
the genome. The score is minus the change in mean log-likelihood (variant minus reference), so a higher score means
a more damaging variant. Scoring runs on one H100 with the first run's container image, pinned by digest, and its
scoring code (`../rare-disease/evo2/launch.sh`, `score_windows.py`). The checkpoint is `arcinstitute/evo2_7b` on
Hugging Face at revision bda0089f (config `evo2-7b-1m.yml`). The image is
`ghcr.io/rewire-bio/evo2-cloud-operations@sha256:22970720ff8b...`, with scoring code from evo2-cloud-operations
commit f8ae11c. Weights are in bfloat16, with FP8 input projections and flash attention, and TF32 off
(`evo2/run2/receipt.json`). Variants the first run had already scored, with the same checkpoint, image and code, are
reused, from its `run/scores.tsv`, `run/scores_extra.tsv` and `run3/scores.tsv`, taking the first file that has a
key. A random 200 of them are scored again in the new run to check that
the reused scores reproduce; where both exist, the reused score is used. The 200 are
`random.Random(20261111).sample(keys, 200)` over the first run's 8,592 scored keys sorted as text, and are listed in
`evo2/check_keys.tsv`. Evo 2 was added after
every other result had been seen (section 2.11), so it is exploratory.

## Baselines

**Lookups combined.** A logistic regression on the four lookups, the model class, and score-by-class terms, as in
the first run. Scores are standardised (scikit-learn `StandardScaler`) and the model is
`LogisticRegression(max_iter=2000, class_weight="balanced")` [@scikitlearn2026linear]. A missing score is set to 0
on its raw scale before standardising, and an indicator for the missing score is added (SNVs only, and only when
some SNV lacks that score). Class indicators use "other" as the reference, and each score has an
interaction term with each class. The label is whether a variant is the diagnostic. Within a case, variants are
ranked by predicted probability. "Lookups + X" adds one score X, with its own class interactions, to the same
features. Adding a score never changes another model. In the indel study the only lookup is our SpliceAI score.

**Rule-based.** Each variant gets ACMG/AMP evidence points [@richards2015standards] from the rule or tool for its
class, and variants are sorted by points (Table 5). This is what a lab applying ClinGen calibrations would do with
the same tools.

| Variants | Points |
|---|---|
| Nonsense SNVs, frameshift indels, canonical splice SNVs and indels | +8 (PVS1 tier) |
| Missense SNVs | the higher of the REVEL level and the SpliceAI level |
| Inframe indels | the higher of +2 (PM4) and the SpliceAI level |
| Everything else | the SpliceAI level |

Table: Points in the rule-based baseline.

The REVEL level follows Bergquist et al. [@bergquist2025calibration], their Table 1: +4 at 0.932 or more, +3 at 0.879, +2
at 0.773, +1 at 0.644; -1 at 0.290 or less, -2 at 0.183, -3 at 0.052, -4 at 0.016; otherwise 0. These match Pejaver
et al. [@pejaver2022calibration], whose very strong benign level (0.003 or less) is merged into -4 here. Where REVEL
is missing, the AlphaMissense level from the same table is used: +4 at 0.990, +3 at 0.972, +2 at 0.906, +1 at 0.792;
-1 at 0.169 or less, -2 at 0.099, -3 at 0.070; otherwise 0. ClinGen asks for one tool chosen in advance, so REVEL is
the tool and AlphaMissense fills gaps. The SpliceAI level follows Walker et al. [@walker2023using], their Table 1: +1 at a
maximum delta of 0.2 or more, -1 at 0.1 or less, 0 between. Walker calibrates these to moderate but recommends
applying them at supporting weight, as here. CADD has no calibration outside missense, so it sets no level. A variant
whose class needs a score and has none gets no level and ranks below every scored variant.

**Ties.** The rule is not a total order, because scores from different tools are not on one scale. Each SNV level has a
tool: for missense, REVEL (or AlphaMissense where REVEL is missing), or SpliceAI when the SpliceAI level is strictly
higher or neither REVEL nor AlphaMissense has a score; SpliceAI for the other non-PVS1 classes; none for the PVS1
tier. A background SNV with the same points as the diagnostic competes with it if both levels came from the same tool and the background's score is higher, or the scores are equal and its
CADD PHRED is at least as high. Otherwise (different tools, or no tool) it competes if its CADD PHRED is at least as
high. Indels at the same points are compared on our SpliceAI score, inframe indels included; equal or missing on both
counts against the diagnostic. Competitors are counted by comparing each background variant with the diagnostic.

**Rule + AVI.** The same points, then AVI within a points level (a missing AVI ranks below a scored one), then the tie
rule above.

The rule does not check whether loss of function is the disease mechanism for a gene, or whether a nonsense or
frameshift variant escapes nonsense-mediated decay.

## Models and test sets

The SNV (200) and indel (120) studies are analysed with 5-fold cross-validation grouped by person (scikit-learn
`GroupKFold`), so no case is scored by a model that saw its person. These results are descriptive.

The replication study is the test set. Every fitted model is fitted once on all 200 SNV cases under the main filters
(grpmax below 0.1% and the inheritance filter) and applied unchanged to the 500 ("locked"). The two exceptions are
the 1% filter check and the EVEE analysis, whose locked models are fitted on the 200 under the 1% filter. The rule
scores have nothing to fit. As a descriptive check, the models are also cross-validated on the 500.

The EVEE analysis covers the 200 and 500 SNV cases. It uses the 1% filter and restricts each list to variants with
an EVEE row; cases whose diagnostic has none are dropped. EVEE, every single score and Lookups combined (fitted on
the full 1% lists) are ranked on those lists. Its result is labelled as agreement with
ClinVar labels on in-sample variants, not diagnostic performance.

## Outcomes and the reading rule

For each score and case, the competitors are the background variants that score at least as high as the diagnostic
(ties count against the score; for the rule scores, the pairwise comparison above). Top *k* is the share of cases with
fewer than *k* competitors. Top 1 is the share of cases where the diagnostic ranks first, and is the outcome for
every pre-specified contrast. Gains are in percentage points.

The primary contrast is top 1 of Lookups + AVI minus top 1 of Lookups combined, both locked, on the 500 replication
cases, under the main filters, with a paired 95% interval. The two models differ only by AVI. The key secondary
contrast is Rule + AVI minus Rule-based on the same cases: does AVI help the stronger per-type baseline, or only make
up for the fitted model's weakness?

The smallest gain of interest is 5 points at top 1: a smaller gain would not change what a lab should do. The reading
rule, fixed before the replication cases were scored:

- "AVI adds a gain that matters" if the primary interval lies above 0, the estimate is at least 5 points, and "the
  key secondary points the same way" (read here as its estimate above 0).
- "No gain that matters" if the primary interval lies inside -5 to +5.
- Otherwise "inconclusive".

The secondary family is GPN-MSA, GPN-Star V and GPN-Star M, each added to the lookups combined, locked, on the 500. A
secondary score is said to add something only if its Bonferroni-adjusted 98.33% paired interval excludes zero.

## Statistics

Intervals are percentile intervals from 1,000 bootstrap resamples of people, each person drawn with all their cases.
Each stratum (all cases, each diagnostic group, each consequence class within the non-coding group, and the gnomAD
found and absent strata) has its own draws. Each draw is `rng.choice(n, size=(1000, n), replace=True)` over the
stratum's people sorted by ID, with the generator seeded as follows:

```
seed = int.from_bytes(sha256(f"{study}/{analysis}/{stratum}").digest()[:8],
                      "big")
rng = np.random.default_rng(seed)
```

The draws are made before any score is evaluated and reused by every score and every
paired difference in that stratum, so adding or dropping a score never moves another interval. A rate in a draw is
the number of hits over the number of cases across the drawn people.

The two checks that compare subsets use one set of draws over all people instead: the gnomAD found and absent split
(label `<study>/<analysis>/gnomad split`) and the training-chromosome split (`<study>/heldout/splits`). Both subsets
and their difference are computed on those draws, so the difference is paired by person; a draw with no case in a
subset is left out of that interval. The found and absent strata above are used only for single-score top *k*.

For locked models the intervals are conditional on the models fitted on the 200; they do not include the variation
from refitting on a different training set. Cross-validated intervals are conditional on one fold assignment.

Every gain has a paired interval: the per-case difference in the top-*k* indicator, resampled with the same draws.
Each paired gain also has a two-sided bootstrap p-value against zero: twice the smaller share of draws at or below
zero and at or above zero, floored at 1/1,000 and capped at 1. This p-value is the percentile interval inverted,
not a test under a null distribution. With its floor at 0.001, no descriptive result can be checked against a
correction for the number of contrasts. The p-values were added after the results without
Evo 2 and change no interval. They are unadjusted. Outside the primary contrast and the GPN family they are
descriptive, like the intervals. A p-value of 0.001 is the floor and means 0.001 or less. With 1,000 draws, tail
intervals beyond 95% rest on few draws.

The GPN family has Bonferroni-adjusted 98.33% intervals (100 minus 5/3). The Evo 2 contrast has a 95% interval and,
for comparison with the GPN family, a 98.75% interval (the Bonferroni level for four scores). The GPN rows keep their
pre-specified 98.33%, so the four together are not held to a 5% family-wise error.

Everything else is descriptive, with unadjusted 95% intervals: top 5 and top 10; each diagnostic group with at least
20 cases, and the non-coding group by CADD class (groups under 20 cases are marked small); Lookups combined minus
Rule-based and Lookups + AVI minus Rule-based; the cross-validated 200 and 500; the indel study; the checks in
section 2.10; each single score's top 1 by class; EVEE; and Evo 2.

Analyses ran with Python 3.12.14, NumPy 2.5.4, pandas 3.1.0, scikit-learn 1.9.1 and SciPy 1.18.1, pinned in
`requirements.txt`.

## Checks

Each check reruns the analysis with one change and the locked models unchanged, except where stated.

- AVI training chromosomes. AVI was trained on gnomAD frequency classes on chromosomes 1, 4, 7, 8, 10, 13 and
  15, from a gnomAD that includes the 1000 Genomes people. Cases are split by whether the diagnostic is on one of
  those chromosomes, and each case is ranked only among its background variants on chromosomes of the same split.
  The models are fitted on the full lists. The primary contrast and the key secondary are reported in each split,
  with the difference between splits. If the gain is larger on training chromosomes and the interval of the
  difference excludes zero, the primary result is reported as possibly inflated by training overlap.
- Looser background filter. grpmax below 1%, with locked models refitted on the 200 under the same filter. This
  tests how much a denser background lowers top 1. It also matches the 1% threshold Genomics England applies to
  recessive modes (section 2.3), though it applies 1% to every gene.
- No background inheritance filter.
- All of 1000 Genomes left out. Every released 1000 Genomes sample is subtracted instead of the family, using the
  adj counts of all 2,500. Where the subset file has no row, each group's AN is reduced in proportion to that group's
  share of 1000 Genomes samples among its gnomAD genomes. This is the literal fix proposed in review; it also removes
  people of the same population who would be in gnomAD for a real patient, and about a fifth of gnomAD's eas and sas
  genomes.
- Diagnostic in gnomAD or not. In the first run's 200 cases, 76% of diagnostics were absent from gnomAD, while
  leaving the family out leaves about 2% of all background variants absent, counted before the frequency filter
  (`results/gnomad_overlap.md`). In this run,
  386 of the 500 diagnostics are absent (section 3.1). A score that learned
  gnomAD presence could pick the diagnostic for that reason alone. The primary and key secondary contrasts are
  reported separately for cases whose diagnostic is in full gnomAD (AC above 0) and cases whose diagnostic is absent,
  with the difference (absent minus found) on the shared draws described in section 2.9. If a gain is larger in the absent cases and that interval
  excludes zero, it is reported as possibly due to gnomAD presence.
- Ties. The share of cases in which a background variant ties the diagnostic under the rule.
- Invariance. `tests/test_invariance.py` reruns an analysis with one added score dropped and checks that no row
  not involving it changes.

## Deviations from the protocol

The protocol was revised after a design review and a code review. Table 6 lists the changes that alter the design
or the numbers. All are dated 2 November 2026. The design-review changes were made before any result of this run.
The code-review changes were made after the cross-validated results on the 200 and 120 had been seen and before the
500 were scored. The last two were made after the results without Evo 2 had been seen. The protocol planned a tagged
freeze with case-file checksums before the 500 were built; it was dropped, so the order of events rests on the git
history (commits 6624e6a, d14fc58, b71d33b).

| Section | Change | Reason |
|---|---|---|
| 2.2 | Leave out the person and family, not all of 1000 Genomes; the latter becomes a check | Neither gnomAD release carries subtractable 1000 Genomes counts; leaving out the family matches a real patient |
| 2.2 | Haploid adj DP is 10, not 5 | Reproduces 1,592 of 1,592 chrX rows outside the PARs; DP 5 reproduces 4 |
| 2.2 | Allele balance is alternate AD over DP | Reproduces 662 of 662 chr22 rows, against 656 for alternate over summed AD |
| 2.6, 2.8 | Rule-based baseline uses ACMG points; primary contrast is Lookups + AVI against Lookups combined | The first draft's primary mixed the model-versus-rule difference with AVI's gain |
| 2.1 | 50 new replication people | The first run's replication people and results had already been seen |
| 2.3 | Inheritance mode added; GNAS (unknown) diagnostics heterozygous | About half of the first run's diagnostics were carrier states |
| 2.3 | The two-variant count includes the person's variants of the other type | A pipeline sees SNVs and indels together; counting one type dropped 32 background SNVs in 28 of the 200 cases |
| 2.1 | Draw walks a shuffled list with the gene cap; MSH3 and MUTYH excluded | A fixed-size sample cannot apply the cap; those genes cause cancer predisposition |
| 2.1 | Background P/LP removal covers every review status | 0 and 1 star P/LP variants left a second plausible answer in 12 SNV and 7 indel cases |
| 2.1 | Indels with REF and ALT both longer than 1 bp left out | 9 of 120 indel diagnostics were such; backgrounds hold none; SpliceAI cannot score them |
| 2.6 | Rule ties written out: different tools never compared, CADD between them | Comparing REVEL with SpliceAI values moved rule-based top 1 by 4 points on the 200 |
| 2.7 | Locked models fitted under the main filters, except the 1% check and EVEE | The code had refitted under each check's filters |
| 2.10 | Check added: gains split by whether the diagnostic is in gnomAD | In the first run's cases, leaving the family out leaves 2% of background variants absent against 76% of diagnostics |
| 2.11 | Planned freeze (tag and case-file checksums before the 500 were built) dropped | The order of events is recorded in the git history instead |
| 2.5, 2.9 | Evo 2 scored on all cases, exploratory, 98.75% interval (after the results without Evo 2) | The author's decision after tag `amp-2026-v5` |
| 2.9 | Two-sided bootstrap p-values added (after the results without Evo 2) | To state each finding at its significance level; no interval changes |

Table: Deviation log entries that change the design or the numbers. The full log, with dates, is section 13 of
`PROTOCOL.md`.

The protocol also drops the first run's inherited cancer panel and keeps GPN-Star as a secondary score, both
decisions of the author.

## Data and code availability

Code, case files, per-case scores (except those whose licences forbid redistribution), result files and logs are in
the repository `github.com/rewire-bio/amp-2026-studies` (MIT licence), directory `rare-disease-v2/`
[@rewire2026rdv2]. Tag `amp-2026-v5` holds every result without Evo 2; tag `amp-2026-v6` holds the results
with Evo 2 and this paper. Every number in the results comes from `results/*.json` and `results/*.md` at those tags.
P-values exist only at `amp-2026-v6`, from the rerun that adds Evo 2; every other number is unchanged from
`amp-2026-v5`. Figure 3 shows the scripts in the order they run.

![The scripts, in the order they run.](figures/diagram_pipeline.pdf){width=80%}

CADD and REVEL are free for non-commercial use but not licensed for redistribution, so their per-case scores are not
committed; `lookup_scores.py` fetches them again. EVEE's table is CC BY 4.0. The other inputs are public: PanelApp,
ClinVar, MANE, the 1000 Genomes high-coverage panel, gnomAD v3.1.2 and its HGDP and 1000 Genomes subset,
AlphaMissense, the AlphaGenome AVI table, and the GPN-MSA and GPN-Star tables. The `README.md` in the study
directory gives the order in which to run the scripts.

```{=latex}
\clearpage
```

# Results

Table 7 lists the main findings. All are top-1 gains in percentage points on the 500 replication cases with locked
models, unless the row says otherwise. Each interval is at the level set for its row in section 2.9. A row is
significant if its interval excludes zero at that level; a bound of exactly 0.0 does not exclude it. The last column
applies the 5-point rule of section 2.8 to each row: Yes if the row is significant and the estimate is at least 5
points, No if the interval lies inside -5 to +5, otherwise Inconclusive. Only the primary row has a verdict under the
reading rule. For every other row the column describes the estimate; it is not a finding, and the descriptive rows
are unadjusted.

| Finding | Estimate | Interval | Level | Significant | Matters |
|----------------------------|----------|--------------|--------------------|------------|--------------|
| Q1: AVI added to the lookups (primary) | +1.8 | -0.4 to +4.0 | 95%, primary | No | No (verdict) |
| Q3: AVI added to the rule (key secondary) | +2.6 | +0.4 to +5.6 | 95%, unadjusted; no test pre-specified | Excludes zero; not a pre-specified test | Inconclusive |
| Q2: lookups minus rule | +3.4 | 0.0 to +7.4 | 95%, descriptive | No | Inconclusive |
| Q2: lookups minus rule, non-coding | +12.0 | +4.7 to +20.9 | 95%, descriptive | Yes | Yes, in this group |
| Q1: AVI added to the lookups, non-coding | +5.1 | +0.6 to +10.5 | 95%, descriptive | Yes | Yes, in this group |
| Q3: AVI added to the rule, non-coding | +8.6 | +3.2 to +15.1 | 95%, descriptive | Yes | Yes, in this group |
| Lookups + AVI minus rule (mixes AVI with the model-versus-rule difference) | +5.2 | +1.6 to +9.2 | 95%, descriptive | Yes | Not applied |
| Q1 on the 200 SNV cases, cross-validated | +4.0 | -2.0 to +9.5 | 95%, descriptive | No | Inconclusive |
| Q4: GPN-MSA added to the lookups | -0.2 | -1.6 to +1.4 | 98.33%, secondary | No | No |
| Q4: GPN-Star V added to the lookups | +0.4 | -1.8 to +2.4 | 98.33%, secondary | No | No |
| Q4: GPN-Star M added to the lookups | -0.2 | -1.6 to +1.2 | 98.33%, secondary | No | No |
| Evo 2 added to the lookups | +0.8 | -0.6 to +2.4 | 98.75%, exploratory | No | No |
| Evo 2 added, 200 SNV cases, cross-validated | +1.0 | -4.5 to +6.0 | 98.75%, exploratory | No | Inconclusive |
| Evo 2 added, 120 indel cases, cross-validated | +1.7 | -9.0 to +12.8 | 98.75%, exploratory | No | Inconclusive |
| Evo 2 added, 120 indel cases, cross-validated, top 5 | +4.2 | +0.8 to +9.5 | 98.75%, exploratory | Yes | Inconclusive |

Table: Main findings. Lookups: the lookups combined by the fitted model. Rule: the rule-based ranking. Non-coding:
the 175 non-coding diagnostics, one of four diagnostic groups. The row that mixes two changes is the contrast the
deviation log removed from the primary (Table 6); the 5-point rule is not applied to it. The six descriptive rows are
among the 944 paired top-1 contrasts in the 16 result files (808 without Evo 2), each descriptive one at an
unadjusted 95%. 175 of the 944 have p below 0.05; if there were no real differences, at most about 47 would, by
chance.

## Cases

Every quota was filled in all three studies, and every build check passed: each diagnostic has grpmax below 0.1%,
each background variant below 1%, no gene exceeds 10 diagnostics in a group, and each position is used once
(`results/build_snv.md`, `results/build_indel.md`, `results/build_replication.md`). The replication draw
excluded 852 positions used as diagnostics before. Under the main filters no case lost its diagnostic.

**People.** The first run's 50 people (24 male, 26 female) are afr 12, nfe 11, eas 11, sas 10 and amr 6 by gnomAD's
population inference, which the subset's sample metadata gives for all 50, including the 16 not in the release. The 50 replication people (26 male, 24 female) are in afr 10, eas 10, nfe
10, sas 10, amr 9 and fin 1; all are in the release (`data/people.tsv`).

**Diagnostic mix.** Table 8 gives the groups. The replication diagnostics span 256 genes. The gene cap was reached
only in the non-coding groups (10 diagnostics from one gene in the 200 and in the 500); elsewhere the largest gene
gave 2 to 8. The non-coding SNV pool was the tightest: 214 positions were still available for the 500's 175.

| Study | Group | Diagnostics | Genes | Largest gene |
|---------------|---------------|---------------------------------|------|--------|
| SNV (200) | missense | 60 | 50 | 4 |
| | canonical splice | 40 | 33 | 4 |
| | nonsense | 30 | 26 | 3 |
| | non-coding | 70 (intronic 38, splice region 23, synonymous 6, UTR 3) | 46 | 10 |
| Indel (120) | frameshift | 40 | 34 | 3 |
| | inframe | 30 | 26 | 3 |
| | canonical splice | 30 | 27 | 2 |
| | non-coding | 20 (intronic 13, UTR 7) | 14 | 5 |
| Replication (500) | missense | 150 | 110 | 6 |
| | canonical splice | 100 | 76 | 8 |
| | nonsense | 75 | 60 | 5 |
| | non-coding | 175 (intronic 84, splice region 63, synonymous 15, UTR 13) | 117 | 10 |

Table: Diagnostic variants by group. Largest gene is the most diagnostics from one gene in the group (cap 10).

About half or more of the diagnostics are homozygous (Table 9): 105 of the 200, 74 of the 120 and 293 of the 500,
nearly all in biallelic genes. Diagnostics in full gnomAD (AC above 0) are a minority: 114 of the 500 (386 absent), 52 of the 200 (148 absent)
and 29 of the 120 (91 absent).

| Study | Monoallelic | Biallelic | X-linked | Unknown | Het | Hom | Hemi |
|------------------|----------|---------|--------|--------|-----|-----|-----|
| SNV (200) | 71 | 96 | 33 | 0 | 77 | 105 | 18 |
| Indel (120) | 43 | 73 | 4 | 0 | 44 | 74 | 2 |
| Replication (500) | 177 | 277 | 44 | 2 | 187 | 293 | 20 |

Table: Diagnostics by mode of inheritance of their gene and by genotype. X-linked combines monoallelic in females
(SNV 14, indel 2, replication 16) and biallelic in females (19, 2, 28). Monoallelic includes PanelApp's "both".

**Lists.** The two filters shorten the lists (Table 10). On the 500, 80,870 case rows (SNVs and indels) fall to
20,580 after the frequency filter and 11,856 after the inheritance filter; dropping the 1,381 background indels that
counted only toward the two-variant rule leaves 10,475 SNVs ranked, a median of 19 per case with the diagnostic.
Background variants removed as ClinVar pathogenic or likely pathogenic at build were 52 (200), 32 (120) and 110
(500).

| Study | Case rows | After frequency | After inheritance | Ranked | Median list |
|------------------|--------|---------|----------|--------|--------|
| SNV (200) | 36,508 | 10,956 | 6,189 | 5,401 | 26 |
| Indel (120) | 21,042 | 6,410 | 3,615 | 576 | 4 |
| Replication (500) | 80,870 | 20,580 | 11,856 | 10,475 | 19 |

Table: Variants under the main filters. Case rows hold SNVs and indels; Ranked keeps the study's own type. The
median list includes the diagnostic.

The checks change list length as intended (Table 11). Without the inheritance filter the median list on the 500 is
34.5 with the diagnostic; the first run's grpmax side check left a median list of 34.5 variants, diagnostic included,
on the 498 of its 500 cases that kept their diagnostic.
The looser 1% filter gives 69.5. Leaving out all of 1000 Genomes instead of the family adds little (21).

| Analysis | SNV (200) | Indel (120) | Replication (500) |
|---|---|---|---|
| Main | 26 | 4 | 19 |
| Background grpmax below 1% | 78.5 | 14 | 69.5 |
| No background inheritance filter | 51.5 | 7 | 34.5 |
| All of 1000 Genomes left out | 25.5 | 4 | 21 |

Table: Median list length per case, diagnostic included, by analysis.

**Coverage.** CADD, AVI, the three GPN scores and Evo 2 cover every SNV in the lists (Table 12). SpliceAI covers 83% of
ranked SNVs (8,729 of 10,475 on the 500) and 484 of the 500 diagnostics. REVEL scores 168 diagnostics and
AlphaMissense 143, against 150 missense diagnostics. Both are looked up over all transcripts, while the class comes
from the MANE transcript, so a variant can have a REVEL score without being missense, and a missense diagnostic can
lack an AlphaMissense score. Of the lookups, only our SpliceAI run applies in the indel study; it scores 483 of 576 ranked indels and
111 of the 120 diagnostics, so 9 diagnostics rank last for SpliceAI.

| Score | SNV (200): variants | diagnostics | Replication (500): variants | diagnostics |
|---|---|---|---|---|
| CADD | 5,401 | 200 | 10,475 | 500 |
| SpliceAI | 4,485 | 192 | 8,729 | 484 |
| REVEL | 1,244 | 68 | 2,090 | 168 |
| AlphaMissense | 1,025 | 57 | 1,885 | 143 |
| AVI | 5,401 | 200 | 10,475 | 500 |
| GPN-MSA, GPN-Star V, GPN-Star M | 5,401 | 200 | 10,475 | 500 |
| Evo 2 | 5,401 | 200 | 10,475 | 500 |

Table: Variants and diagnostics with a score in the ranked lists, main filters. The lists hold 5,401 SNVs (200
cases) and 10,475 (500 cases).

Evo 2 also scores all 576 ranked indels and all 120 indel diagnostics, where SpliceAI scores 483 and 111.

The model classes, set from the 200 SNV cases, merge synonymous (6 diagnostics) and UTR (3) into other; from the 120
indel cases, UTR (7) is merged. No diagnostic's class differs between the build and the analysis, and no case
variant lacks a MANE row in CADD.

## Primary result

Adding AVI to the combined lookups does not raise top 1 by a gain that matters. On the 500 replication cases with
locked models, the lookups combined rank the diagnostic first in 90.6% of cases (87.8 to 93.4) and Lookups + AVI in
92.4% (89.8 to 94.8). The gain is +1.8 points (-0.4 to +4.0), p = 0.162. The interval lies inside -5 to +5, so under
the reading rule the result is "no gain that matters" (Figure 4).

![Top-1 gain from AVI, with paired 95% intervals. The dashed lines mark no gain and plus or minus 5 points, the
smallest gain of interest. AVI scores SNVs only, so the indel study has no AVI contrast.](figures/gains_primary.pdf){width=90%}

Top 1 is the pre-specified outcome. At top 5 (descriptive) the gain is +1.0 (+0.2 to +2.2), p = 0.028; at top 10 it
is +0.2 (-0.4 to +0.8). Without AVI the locked lookups already place the diagnostic in the top 5 in 97.0% of cases
(95.4 to 98.4), so there is little to gain there either.

The other estimates point the same way (95%, descriptive). Cross-validated on the 500, the lookups reach 92.8%
(89.8 to 95.6) and the gain is +0.8 (-1.8 to +3.2). On the 200 SNV cases, cross-validated, the gain is +4.0 (-2.0 to
+9.5), p = 0.178.

By diagnostic group (95%, descriptive, four groups), the gain is in the 175 non-coding diagnostics: +5.1 (+0.6 to
+10.5), p = 0.030, from 84.0% (77.4 to 89.9) to 89.1% (85.0 to 93.3). On the 150 missense diagnostics it is 0.0 (-5.5 to +4.6). On canonical splice
and nonsense diagnostics both models rank every diagnostic first, so the gain is 0.0. Within the non-coding group the
gain is +4.8 (0.0 to +11.6) for intronic and -1.6 (-8.8 to +5.4) for splice region diagnostics, and +30.8 (+8.3 to
+50.0) for the 13 UTR diagnostics, a small group.

## How often the lookups rank the diagnostic first

The lookups alone rank the diagnostic first in about 9 cases in 10. On the 500, Lookups combined (locked) reaches
90.6% (87.8 to 93.4) and the rule-based ranking 87.2% (82.8 to 91.0). The difference is +3.4 (0.0 to +7.4), p =
0.054 (95%, descriptive). At top 5 they reach 97.0% (95.4 to 98.4) and 97.8% (96.6 to 98.8). Combining the lookups
does far better than any one of them: the best single lookup over all cases is CADD at 61.0% (53.4 to 68.8), because
each lookup is strong only in its own variant classes (Table 13, Figure 5).

| 500 cases | Cases | Lookups combined (locked) | Rule-based | Lookups minus rule | Best lookup |
|-------------|-------|------------------|------------------|------------------|-----------------------|
| Canonical splice | 100 | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 0.0 (0.0 to 0.0) | SpliceAI 96.0 (91.1 to 100.0) |
| Nonsense | 75 | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | 0.0 (0.0 to 0.0) | CADD 98.7 (95.5 to 100.0) |
| Missense | 150 | 87.3 (82.1 to 92.5) | 90.0 (84.0 to 95.4) | -2.7 (-7.3 to +2.8) | REVEL 90.7 (84.7 to 95.7) |
| Non-coding | 175 | 84.0 (77.4 to 89.9) | 72.0 (61.7 to 81.1) | +12.0 (+4.7 to +20.9) | SpliceAI 86.3 (81.3 to 91.3) |
| All | 500 | 90.6 (87.8 to 93.4) | 87.2 (82.8 to 91.0) | +3.4 (0.0 to +7.4) | CADD 61.0 (53.4 to 68.8) |

Table: Top 1 by diagnostic group on the 500, with 95% intervals (descriptive; the difference is one of five strata).
Best lookup is the single lookup with the highest top 1 in the group.

Within each group the combined model is now close to the best single lookup, unlike in the first run. The rule falls
behind on non-coding diagnostics, +12.0 (+4.7 to +20.9), p = 0.001, where it has only SpliceAI's supporting level to
go on; it ranks none of the 13 UTR diagnostics first. On missense diagnostics the rule is ahead on the estimate at top
1, -2.7 (-7.3 to +2.8), and at top 5, -2.7 (-5.4 to -0.6), p = 0.024.

```{=latex}
\begin{figure}[p]
\centering
\makebox[\textwidth][c]{\rotatebox{90}{\includegraphics[width=0.92\textheight]{figures/top1_by_group_replication.pdf}}}
\caption{Top 1 by diagnostic group on the 500 replication cases, locked models, with 95\% intervals. Groups under 20
cases are marked small.}
\end{figure}
```

On the 200 SNV cases (cross-validated, descriptive), the lookups reach 84.0% (75.5 to 91.5) and the rule 81.0% (73.5
to 88.0), a difference of +3.0 (-4.5 to +9.5). The rule is ahead on canonical splice diagnostics, 97.5% (91.7 to
100.0) against 82.5% (68.3 to 95.6), a difference of -15.0 (-29.3 to -2.5), p = 0.034. So on the 200 SNV cases the
fitted model still loses to the rule-based ranking in one variant type, as in the first run. The lookups lead on
non-coding diagnostics, +17.1 (+6.2 to +28.3). On the 120 indel cases, ranked among a median of 4 indels, the lookups
(our SpliceAI score with the indel classes) reach 84.2% (75.4 to 92.4) and the rule 76.7% (65.8 to 86.8), a
difference of +7.5 (+3.2 to +12.6), p = 0.001.

## Rule-based ranking plus AVI

On the 500, Rule + AVI minus Rule-based is +2.6 (+0.4 to +5.6), p = 0.034, at an unadjusted 95%. The protocol sets
no test for this contrast and uses it only for its direction, so this is not a confirmed gain. The estimate is under
the 5-point smallest gain of interest and the interval reaches above it, so by the 5-point rule it is inconclusive.
The primary interval decides the verdict.

The gain is in non-coding diagnostics: +8.6 (+3.2 to +15.1), p = 0.001, from 72.0% (61.7 to 81.1) to 80.6% (72.7 to
87.6) (95%, descriptive, one of four groups). Within them it is +13.1 (+3.2 to +25.0) for intronic and +3.2 (0.0 to
+8.8) for splice region diagnostics. AVI orders variants that the rule puts on the same points level, and for non-coding diagnostics that
level is mostly SpliceAI's supporting level. On missense diagnostics the gain is -1.3 (-3.8 to +1.3); canonical splice
and nonsense diagnostics are at 100% either way.

On the 200 SNV cases (cross-validated, descriptive) the gain is +4.5 (+0.5 to +8.0), p = 0.028, and +11.4 (+2.5 to
+21.6) on non-coding diagnostics. In the main analysis of the 200 it sits in cases whose diagnostic is absent from
gnomAD: +6.1 (+1.9 to +11.5), against 0.0 (-6.7 to +5.0) in cases where it is found. The difference, +6.1 (-0.6 to
+13.8), does not exclude zero, so it is not flagged. Under the all-1000-Genomes check it is (section 3.8).

## GPN-MSA and GPN-Star

None of the three GPN scores adds a gain that matters. Added to the lookups, locked, on the 500, with Bonferroni-adjusted
98.33% intervals (Figure 6):

- GPN-MSA: -0.2 (-1.6 to +1.4), p = 0.844.
- GPN-Star V: +0.4 (-1.8 to +2.4), p = 0.666.
- GPN-Star M: -0.2 (-1.6 to +1.2), p = 0.886.

Each interval includes zero and lies inside -5 to +5. Alone, they rank the diagnostic first in 51.2% (43.8 to 58.2),
59.4% (51.8 to 66.0) and 55.8% (48.0 to 63.4) of the 500 cases, against CADD's 61.0% (53.4 to 68.8) and AVI's 83.6% (78.2 to
88.8). On non-coding diagnostics, where AVI alone reaches 71.4% (62.2 to 80.2), they reach 28.6% (19.9 to 37.7),
32.0% (22.6 to 41.7) and 30.9% (20.8 to 41.5). Cross-validated on the 200 (95%, descriptive) the gains are -2.0 (-7.0
to +2.5), -0.5 (-5.5 to +4.0) and +0.5 (-4.0 to +4.5).

In the checks (section 3.8), 1 of the 36 all-case GPN top-1 contrasts excludes zero at 98.33%: GPN-Star M on the 200
with all of 1000 Genomes left out, +5.5 (+1.0 to +10.5), p = 0.001. On the 500 under the same check, locked, it is
0.0 (-1.4 to +1.4).

![Top-1 gain from each added score, with paired 95% intervals unless the label gives another level. Dashed lines mark
no gain and plus or minus 5 points.](figures/gains_added.pdf){width=85%}

## Evo 2 (exploratory)

**The run.** Evo 2 scored every case variant in the three studies: 16,173 distinct SNVs and indels. Of these, 8,592
scores were reused from the first run and 7,581 were new. The new run scored 7,781 variant windows (the 7,581 and the
200 reused variants scored again as a check) and 7,760 reference windows, 15,541 in all, in 5,343 s on one H100 (2.9
windows per second, peak memory 16.8 GiB). The pod ran from 14:00 to 15:42 UTC on 2 November 2026 at USD 3.99 per
hour, about USD 6.80 in all (`evo2/run2/receipt.json`, `evo2/run2/events.jsonl`). The 200 rescored variants (179 SNVs
and 21 indels, `evo2/check_keys.tsv`) gave scores identical to the reused ones: the largest absolute difference is
0.0. Evo 2 covers every ranked variant and every diagnostic (Table 12).

**Results.** On the 500, Evo 2 adds no gain that matters to the lookups (Table 14). Locked, Lookups + Evo 2 minus
Lookups combined is +0.8 (98.75% interval -0.6 to +2.4), p = 0.236. On the 200 SNV cases and the 120 indel cases,
cross-validated, the top-1 results are inconclusive: +1.0 (98.75%: -4.5 to +6.0) and +1.7 (-9.0 to +12.8). On the
indel cases the top-5 gain excludes zero at 98.75% (below); this is exploratory. Alone, Evo 2 ranks the diagnostic first in
58.8% (52.0 to 65.4) of the 500, close to CADD, and in 81.7% (73.7 to 89.2) of the indel cases, above SpliceAI's
58.3% (49.6 to 66.4).

| Study | Evo 2 alone | Lookups combined | Lookups + Evo 2 | Gain, 95% | Gain, 98.75% | p |
|------------------|----------|----------|----------|-----------|-----------|------|
| 500, locked | 58.8 (52.0 to 65.4) | 90.6 (87.8 to 93.4) | 91.4 (88.8 to 94.0) | +0.8 (-0.4 to +2.0) | -0.6 to +2.4 | 0.236 |
| 500, cross-validated | 58.8 (52.0 to 65.4) | 92.8 (89.8 to 95.6) | 93.0 (89.8 to 95.8) | +0.2 (-1.4 to +1.8) | -1.8 to +2.2 | 0.844 |
| 200 SNV, cross-validated | 51.5 (43.5 to 59.5) | 84.0 (75.5 to 91.5) | 85.0 (77.0 to 91.5) | +1.0 (-3.5 to +5.5) | -4.5 to +6.0 | 0.760 |
| 120 indel, cross-validated | 81.7 (73.7 to 89.2) | 84.2 (75.4 to 92.4) | 85.8 (77.7 to 93.5) | +1.7 (-5.9 to +9.9) | -9.0 to +12.8 | 0.746 |
| First run, 120 indel, no filter | 57.5 (47.5 to 67.2) | 61.7 (51.2 to 71.3) | 77.5 (67.2 to 85.5) | +15.8 (+5.0 to +26.8) | | |
| First run, 120 indel, AF below 0.1% | 71.7 (62.4 to 80.5) | 71.7 (61.0 to 80.7) | 80.0 (70.9 to 87.6) | +8.3 (+0.9 to +16.3) | | |

Table: Evo 2, top 1 in percent and gain in points, main filters. Exploratory. The first run's rows come from
`../rare-disease/results/results_indel.json` and `results_indel_gnomad0.001.json`; they have 95% intervals only.

By diagnostic group (95%, descriptive; Figure 7), the largest estimates are +2.9 (0.0 to +6.8) on non-coding
diagnostics of the 500, +5.0 (-5.3 to +17.5) on canonical splice diagnostics of the 200 and +10.0 (-3.9 to +26.7) on
inframe indels; none excludes zero. Within the non-coding group, the 23 splice region diagnostics of the 200 give
+13.0 (-4.3 to +31.8).

At top 5 Evo 2 adds +4.2 on the indel cases (98.75%: +0.8 to +9.5; 95%: +0.8 to +8.2), p = 0.012. All of it is in
the 20 non-coding indels, +25.0 (95%: +5.6 to +44.5); the other three groups gain 0.0. The 98.75% interval also
excludes zero without the inheritance filter, +5.0 (+0.8 to +10.2), p = 0.008, and with all of 1000 Genomes left
out, +4.2 (+0.8 to +9.3), p = 0.010, but not with the 1% filter, +5.0 (0.0 to +11.6), p = 0.024. Each is under 5
points, and top 5 is descriptive. Evo 2 scores every indel, while SpliceAI has no score for 9 of the 120 diagnostics,
which rank last for SpliceAI, and the indel model has no missing-score indicator. Some of the gain may come from
those diagnostics.

Across the four analyses there are 16 all-case Evo 2 top-1 contrasts (the 200 and the 120 cross-validated, and the
500 locked and cross-validated, in each analysis). One has a 98.75% interval that excludes zero: the 200 SNV cases
without the inheritance filter, +5.5 (+0.1 to +12.0). Its p-value, 0.014, is just above 0.0125. The two disagree
because the p-value counts the draws at exactly zero, while the interval's bound is interpolated between a draw at
zero and the next one above it, so this is a borderline result. With 16 contrasts at 98.75%, about 0.2 would exclude
zero by chance. Its top 5 moves the other way, -2.0 (95%: -4.0 to -0.5; 98.75%: -4.5 to 0.0), p = 0.040, which does
not exclude zero at 98.75%.

![Evo 2 added to the lookups, by diagnostic group, with the first run's indel results. The all-case rows have 98.75%
intervals, the others 95%.](figures/evo2_gains.pdf){width=85%}

**Against the first run.** The first run found that Evo 2 added +15.8 points (+5.0 to +26.8) to the indel lookups,
and +8.3 (+0.9 to +16.3) with a gnomAD AF filter below 0.1%. This run finds +1.7 (95%: -5.9 to +9.9). The runs
differ in cases and design, the intervals overlap, and the difference was not tested. Four things changed. In the
first run CADD's precomputed indel table scored 36 of the 120 indel diagnostics but 5,395 of the 5,540 indels in its cases,
because it holds only indels seen in gnomAD; the other 84 diagnostics scored 0 in the lookups,
with no missing-score indicator for indels, while Evo 2 scored all 120. Part of the first run's +15.8 was therefore
Evo 2 filling a gap that the design created. The lists are shorter: the first run ranked a median of 37 indels per
case without a filter and 15 with it; this run ranks 4. The baseline is stronger: the first run's indel lookups were
CADD and SpliceAI and reached 61.7% (51.2 to 71.3); this run drops CADD for indels (section 2.5), adds the indel
classes to the model, and its lookups reach 84.2% (75.4 to 92.4). And the indel diagnostics are a new draw from a
pool without multi-base substitutions.

Within this run the estimates are higher in every check than in the main analysis: +5.8 (98.75%: -0.9 to +13.8)
without the inheritance filter (median list 7), +7.5 (-2.4 to +18.1) with the 1% filter (median 14) and +5.0 (-4.1
to +14.2) with all of 1000 Genomes left out (median 4, as in the main analysis). All four use the same 120 cases,
every 98.75% interval includes zero, and the all-1000-Genomes check gives +5.0 on lists as short as the main
analysis's, so these results do not show a trend with list length.

## EVEE (descriptive)

EVEE ranks the diagnostic first about as often as the lookups, on lists where that comparison says little. With the
1% filter and each list restricted to variants with an EVEE row, 492 of the 500 cases keep their diagnostic (8
diagnostics have no row), with a median list of 26. EVEE alone ranks the diagnostic first in 87.0% (82.5 to 90.9),
the locked lookups in 87.0% (82.2 to 91.0), the cross-validated lookups in 88.0% (82.8 to 92.0) and AVI in 76.0%
(69.7 to 82.0). On the 200 (all keep their diagnostic, median list 29.5), EVEE reaches 81.5% (75.5 to 88.5) and the
lookups 82.0% (74.5 to 88.5). These are single-score 95% intervals; the EVEE files hold no paired contrasts. This is
agreement with ClinVar labels on in-sample variants, not diagnostic performance: EVEE covers ClinVar variants only,
and almost every diagnostic was pathogenic or likely pathogenic in its training labels (section 2.5).

## Checks

**Filters.** Each check changes list length as intended, and AVI's gain to the lookups excludes zero in none of them
(Table 15, 95%, descriptive: three checks, two contrasts each). The denser 1% background lowers top 1 for every
ranking; there the gain from AVI to the rule is +5.2 (+1.4 to +9.6), p = 0.002. Without the inheritance filter, the
gain to the lookups is +4.0 (0.0 to +8.0), p = 0.064. Leaving out all of 1000 Genomes instead of the family changes
little.

| Check (500, locked) | Median list | Lookups combined | AVI added to the lookups | AVI added to the rule |
|------------------------------|--------|------------------|------------------|------------------|
| Main | 19 | 90.6 (87.8 to 93.4) | +1.8 (-0.4 to +4.0) | +2.6 (+0.4 to +5.6) |
| Background grpmax below 1% | 69.5 | 82.0 (76.4 to 87.2) | -1.2 (-6.6 to +3.2) | +5.2 (+1.4 to +9.6) |
| No background inheritance filter | 34.5 | 79.6 (73.2 to 85.4) | +4.0 (0.0 to +8.0) | +2.2 (-1.2 to +5.4) |
| All of 1000 Genomes left out | 21 | 87.8 (83.2 to 91.8) | +0.4 (-2.8 to +3.0) | +2.6 (+0.4 to +5.4) |

Table: Checks on the 500 with locked models. Top 1 in percent and gains in points, 95% intervals. The 1% check uses
models fitted on the 200 under the 1% filter.

**AVI training chromosomes.** 124 of the 500 diagnostics are on the chromosomes AVI was trained on and 376 are not.
Each case is ranked only among background variants on chromosomes of its split, so lists are shorter (median 5 and
13) and top 1 is higher than in the main analysis. The primary gain is +1.6 (0.0 to +4.2) on training chromosomes and
+2.1 (-0.3 to +4.6) on the others. Training minus other is -0.5 (-3.7 to +3.1), p = 0.740, so the primary result is not
flagged as inflated by training overlap. Within this split the lists are short (median 5 and 13) and the locked
lookups reach 93.5% (88.1 to 98.2) on training chromosomes, so the check has little room to show inflation. For the key secondary the gains are +0.8 (-2.6 to +4.6) and +2.4 (+0.5 to +4.7),
and the difference is -1.6 (-5.6 to +2.5). On the 200 (38 and 162 cases), training minus other is -0.4 (-9.8 to +8.6)
for the primary contrast.

**Diagnostic in gnomAD or not.** These results come from the shared draws of the gnomAD split (section 2.9). 114 of
the 500 diagnostics are in full gnomAD and 386 are absent. The primary gain is -0.9 (-6.6 to +5.0) in found cases and
+2.6 (+0.5 to +4.8) in absent cases; absent minus found is +3.5 (-2.2 to +8.8), p = 0.234, so it is not flagged. The point estimates fit a gain that
depends on gnomAD absence, and the check is too imprecise to rule this out. The
key secondary is +3.5 (0.0 to +8.2) in found and +2.3 (0.0 to +5.2) in absent cases, a difference of -1.2 (-5.3 to
+2.8). No difference is flagged in any replication analysis (Figure 8). One flag fires on the 200: under the
all-1000-Genomes check, Rule + AVI minus Rule-based is +6.8 (+2.6 to +12.3) in absent cases and -1.9 (-6.7 to 0.0) in
found ones, a difference of +8.7 (+3.2 to +15.5), p = 0.001. So on the 200 the key secondary gain is possibly due to
gnomAD presence. The 500 do not show it.

![Gain from AVI by whether the diagnostic is in full gnomAD, with paired 95% intervals on shared draws.](figures/gains_gnomad.pdf){width=85%}

**Ties.** Under the main filters no SNV diagnostic ties a background variant under the rule or under Rule + AVI, in
the 200 or the 500. In the indel study 5 of 120 cases have a tie (4.2%): 2 of 40 frameshift and 3 of 30 inframe
diagnostics. With the 1% filter ties occur in 2 of the 500 cases (0.4%) and 11 of the 120 (9.2%).

**Invariance.** `tests/test_invariance.py` reruns seven analyses, each with one added score dropped, and compares
9,571 result rows that do not involve that score: none changed (`logs/cb_invariance.log`). This log is from commit
d14fc58, before Evo 2 and the p-values were added. For those, the check is the comparison with tag `amp-2026-v5`
below. Every value in the 16
result files that does not involve Evo 2 or a p-value is identical to tag `amp-2026-v5`. Adding Evo 2 and the
p-values changed no other number.

**Evo 2 reproducibility.** The 200 reused Evo 2 scores that were scored again in the new run are identical to the
first run's: the largest absolute difference is 0.0 (section 3.6).

# Discussion

For a lab that ranks rare variants with CADD, SpliceAI, REVEL and AlphaMissense after a grpmax filter and an
inheritance filter, adding AVI does not change how often the diagnostic variant comes first by an amount that
matters: +1.8 points (-0.4 to +4.0) on 500 new cases. The lookups already rank it first in about 9 cases in 10, and
in every canonical splice and nonsense case, so there is little room. Where there is room, in non-coding diagnostics,
AVI's estimated gain was larger: +5.1 (+0.6 to +10.5) added to the lookups, locked, and +2.9 (+0.6 to +5.2)
cross-validated; +8.6 (+3.2 to +15.1) added to the rule. These are descriptive results in one of four groups. In this
benchmark, ordering variants within a points level by AVI helped intronic diagnostics (+13.1, +3.2 to +25.0, 95%,
descriptive) but not clearly splice region ones (+3.2, 0.0 to +8.8). This is a subgroup result and would need its
own test before a lab relied on it. GPN-MSA and GPN-Star add no gain that matters, and neither does Evo 2 on SNVs.
On indels the Evo 2 result is inconclusive at top 1, and its top-5 gain, though it excludes zero, is under 5 points
and exploratory.

The rule-based ranking does nearly as well as the fitted model: 87.2% (82.8 to 91.0) against 90.6% (87.8 to 93.4),
a difference of +3.4 (0.0 to
+7.4). It needs no training data and follows published calibrations. Its weak point is non-coding variants, where the
fitted lookups lead by +12.0 (+4.7 to +20.9).

**Why the first run's estimates were larger.** The first run's headline, +6.2 (0.0 to +13.0) on its 500 with
locked models, had no frequency filter, a median list of 208.5 variants per case (diagnostic included), the person
in gnomAD and carrier states counted as diagnoses; its locked lookups reached 61.0% (53.6 to 67.2). Under its grpmax
side check (median list 34.5, on the 498 cases that kept their diagnostic) the lookups reached 78.3% (71.2 to 84.4)
and the gain was +0.2 (-3.6 to +3.8). This run's estimate, +1.8 (-0.4 to +4.0), is close to that side check and
below the headline, though its interval overlaps the headline's; the two runs differ in cases and design, and the
difference was not tested (Figure 9). The overall AF filter of the first run, +7.4 (+1.2 to +14.2), uses gnomAD's
overall frequency, which keeps variants that are common in one group. Two explanations fit these numbers, though
neither was tested. Long lists hold many variants that are common in some group, which a score trained on gnomAD
frequency classes can place below the diagnostic; a grpmax filter removes them before any score sees them. And when
the lookups rank the diagnostic first 9 times in 10, an added score can only gain in the remaining tenth. For Evo 2
on indels, shorter lists (4 against 37) and a stronger baseline (84.2% against 61.7%) are possible explanations, not
shown ones; the first run's baseline also lacked a CADD score for 84 of its 120 indel diagnostics (section 3.6).

![AVI added to the lookups, locked, on 500 replication cases, in the first run under three frequency settings and in
this run. Paired 95% intervals. The first run's grpmax setting has 498 cases, because two cases lost their diagnostic
under that filter.](figures/avi_first_vs_second.pdf){width=90%}

## Limits

- **Short lists.** After the filters the median list is 19 SNVs on the 500 and 4 indels on the 120. These lists hold
  only rare variants in one panel's genes; a lab searching an exome or a genome has more candidates. With the 1%
  filter (median 69.5) the lookups fall to 82.0% (76.4 to 87.2) and the gain from AVI to the rule rises to +5.2 (+1.4 to +9.6).
- **Ceiling.** At 90.6% top 1 there is little room for any added score, and none at all in canonical splice and
  nonsense diagnostics, where both baselines reach 100%. The benchmark cannot show a large gain in those groups.
- **Homozygous diagnostics.** 293 of the 500 diagnostics are homozygous, nearly all in biallelic genes. That fits a
  consanguineous family. In an outbred family most recessive diagnoses are two different variants, which this design
  does not model, because scores rank single variants. The two-variant rule ignores phase. The backgrounds come from
  outbred people, so they lack the many rare homozygous variants that runs of homozygosity give a consanguineous
  patient, and that the biallelic filter would keep. Lists for homozygous diagnostics are therefore shorter than a
  consanguineous patient's.
- **Frequency thresholds.** The main filter applies 0.1% grpmax to every gene. Genomics England uses 1% for
  recessive modes [@genomicsengland2026popfreq], so for the 277 diagnostics in biallelic genes this run filters the
  background more strictly than a Genomics England pipeline would, which shortens their lists and favours every
  baseline. The 1% check matches the recessive threshold: there the lookups fall as in the first bullet, and AVI
  added to the lookups gives -1.2 (-6.6 to +3.2) (95%, descriptive).
- **Diagnostic phenotype and inheritance.** Diagnostic variants were not matched to an intellectual disability
  phenotype or to the inheritance ClinVar records for them, so some heterozygous diagnostics in genes with
  monoallelic or both modes may be carrier states for a recessive condition.
- **ClinVar circularity.** The diagnostics are ClinVar pathogenic or likely pathogenic variants. Some may have been classified
  with computational evidence (PP3) from the same lookups, which favours the lookups and the rule. REVEL's training
  set (HGMD) overlaps ClinVar pathogenic missense variants, and the rule's thresholds were calibrated on ClinVar.
  Variants that reach pathogenic or likely pathogenic with two stars are those with strong evidence, so they are
  easier than a typical unsolved case. Removing every ClinVar pathogenic or likely pathogenic background variant (110
  in the 500) also takes out the strongest competitors a real genome would hold.
- **1000 Genomes people are not patients.** The backgrounds come from people not selected for a rare disease, and the
  diagnostic is inserted, not carried. A patient's genome may hold more rare damaging variants than these.
- **gnomAD presence is only partly fixed.** Leaving out the family leaves about 2% of all background variants in the
  first run's 200 cases absent from gnomAD, counted before the frequency filter, against 386 of the 500 diagnostics
  (77%). A score that learned gnomAD presence can still pick the diagnostic for that reason. The split check did not
  flag the primary contrast, but its gain is entirely in absent cases and the check cannot exclude a presence effect.
  It flagged the key secondary on the 200 under the all-1000-Genomes check.
- **Transcript set.** This run takes CADD's most severe consequence over MANE Select transcripts, and the indel classes
  from MANE Select too. Genomics England's rare disease pipeline annotates consequence against all GENCODE Basic transcripts of a
  gene on Ensembl 90, hosted in CellBase [@genomicsengland2025impact]. Consequences in its high-impact list
  (transcript ablation, splice acceptor, splice donor, stop gained, frameshift, stop lost, initiator codon, start
  lost) in any basic transcript qualify a variant for tier 1. Those in its moderate list (transcript amplification,
  inframe insertion or deletion, missense, splice region, incomplete terminal codon), and non-coding transcript exon
  variants in non-coding genes, qualify it for tier 2, subject to the pipeline's other filters. So a variant that is non-coding on MANE Select here can be
  tier 1 or 2 there through another transcript. Initiator codon variants are high impact there but missense here, and
  stop lost variants fall in "other" here. Some groups, model classes and PVS1 tiers would differ under that
  annotation. This run does not measure how many.
- **Fold sensitivity.** The cross-validated results (the 200, the 120 and the descriptive check on the 500) are
  conditional on one fold assignment, and another split could move them. The locked intervals do not include the
  variation from refitting on another training set.
- **Evo 2 and the p-values came after the results.** Evo 2 was added after every other result had been seen, so it is
  exploratory. Its 98.75% interval and the GPN rows' 98.33% do not hold the four added scores to a 5% family-wise
  error. The p-values were added after the results without Evo 2. They are unadjusted, floored at 0.001, and with
  1,000 draws the tails beyond 95% rest on few draws.
- **Other limits.** Relatedness not recorded in the 1000 Genomes pedigree is not handled. The indel lists are very
  short (median 4), so indel results say little about longer lists. Diagnostics are ranked only among variants of
  their own type; a pipeline ranks SNVs and indels in one list. CADD is not used for indels, although a lab can score
  novel indels with it, so the indel baseline is weaker than a lab's.

# Conclusion

With realistic frequencies and an inheritance filter, the established lookups rank the diagnostic variant first in
about 9 cases in 10, and a rule-based ACMG ranking in nearly as many. Adding AlphaGenome AVI gives +1.8 points (-0.4
to +4.0) on 500 new cases: no gain that matters. AVI's gain was concentrated in non-coding diagnostics
(descriptive, one of four groups: +5.1, +0.6 to +10.5, locked; +2.9, +0.6 to +5.2, cross-validated). GPN-MSA and
GPN-Star add no gain that matters, and neither does Evo 2 on SNVs; on indels the Evo 2 result is inconclusive. The
first run's larger estimates came with longer lists, no grpmax filter and a weaker baseline. A
test on patients' genomes, with longer lists, recessive diagnoses made of two variants and diagnostics not taken from
ClinVar, could change this answer.

# References {-}

::: {#refs}
:::
