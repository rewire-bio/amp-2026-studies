# Research proposal: rare-disease case tests, second run

Status: revised 2 November 2026 after a design review. Design in [`PROTOCOL.md`](PROTOCOL.md), work plan in
[`PLAN.md`](PLAN.md). Tracked in issues [#1](https://github.com/rewire-bio/amp-2026-studies/issues/1) and
[#2](https://github.com/rewire-bio/amp-2026-studies/issues/2).

## Background

The first run ([`../rare-disease/`](../rare-disease/), results at tag `amp-2026-v4`) built cases from a real
1000 Genomes person's rare variants in the PanelApp 285 intellectual disability panel, added one ClinVar
pathogenic variant, and asked where each score ranks it. Its main findings
([`../rare-disease/README.md`](../rare-disease/README.md),
[`../rare-disease/REPLICATION.md`](../rare-disease/REPLICATION.md)):

- The established lookups combined (CADD, SpliceAI, REVEL, AlphaMissense, with consequence) do most of the work.
- AlphaGenome AVI was the one added score whose top-1 gain held on 500 new cases: +6.2 points locked
  (0.0 to +13.0), +7.4 (+1.2 to +14.2) after a gnomAD AF filter of 0.1%.
- Under the stricter grpmax filter, AVI's locked gain on the 500 was +0.2 (-3.6 to +3.8).
- GPN-MSA's +6 point gain on the first 200 cases did not replicate.
- On indels, Evo 2 added +15.8 points at top 1, +8.3 after the AF filter.

## Why a second run

A review of the talk found five problems with that design, and a review of this proposal found a sixth.

1. **The background people are in gnomAD** (issue #1). 1000 Genomes is part of gnomAD, so every background SNV has
   been seen there at least once, while most spiked-in ClinVar variants have not:

   | Variants | n | Not in gnomAD | Below 0.01% | 0.01 to 0.1% | 0.1% or more |
   |---|---|---|---|---|---|
   | 200 cases, spiked-in | 200 | 75% | 24% | 2% | 0% |
   | 200 cases, background | 13,145 | 0% | 17% | 26% | 57% |
   | 500 cases, spiked-in | 500 | 76% | 23% | 1% | 0% |
   | 500 cases, background | 11,208 | 0% | 15% | 29% | 56% |

   In about three quarters of cases, "not in gnomAD" alone picks out the diagnostic variant. A real patient's private
   variants are absent from gnomAD and compete with the causal one. Frequency, and any score that learned it, looks
   stronger than it would in practice. AVI was trained on gnomAD frequency classes.
2. **The realistic filter was a side check.** Published pipelines filter per ancestry group (the Genomics England
   tiering and Exomiser both do). That is the grpmax filter, where AVI showed no clear gain.
3. **The combined model is a weak baseline within each variant type.** On the 200 cases at top 1 it lost to the
   matching single tool in three of the four types: canonical splice 82% against SpliceAI 90%, nonsense 60% against
   CADD 90%, missense 43% against REVEL or AlphaMissense 58%. In the fourth (non-coding or synonymous) it tied
   SpliceAI at 60%. Part of an added score's gain may make up for that.
4. **Two consequence annotations.** Diagnostic variants were grouped by ClinVar's MC field; the model's classes came
   from CADD, whose splice class also takes splice-region variants.
5. **Small groups.** With 30 to 70 cases per type, most per-type differences are noise, and not every gain had a
   paired interval.
6. **Many diagnoses were not diagnoses.** 94 of the 200 and 231 of the 500 diagnostic variants were single
   heterozygous variants in genes that PanelApp lists as biallelic only: carrier states. The backgrounds were not
   filtered by inheritance mode either, as a pipeline would.

## Questions

- **Q1 (primary).** With frequencies a real patient's pipeline would see (gnomAD without the person and their
  family), a grpmax filter below 0.1% and an inheritance filter, does adding AVI to the combined lookups raise top 1
  on 500 cases from 50 people the first run never used?
- **Q2.** How often do the lookups alone rank the diagnostic variant first, overall and by variant type, combined by
  a fitted model and by a rule-based ranking that follows ClinGen guidance?
- **Q3.** Does AVI help the rule-based ranking too, or only the fitted model?
- **Q4 (secondary).** Do GPN-MSA, GPN-Star (vertebrate and mammal) add anything under the same design?

## Scope

In:

- PanelApp 285 (intellectual disability) only.
- 200 SNV cases and 120 indel cases from the first run's 50 people; 500 SNV replication cases from 50 new people.
- Scores: CADD v1.7, SpliceAI, REVEL, AlphaMissense (the lookups); AlphaGenome AVI (primary added score); GPN-MSA,
  GPN-Star V and M (secondary); EVEE (descriptive only, because it covers ClinVar variants only and its scores may
  include its training labels; see the protocol).
- Each case a real diagnosis under its gene's inheritance mode; the background filtered by inheritance as a pipeline
  would.
- A rule-based baseline from ClinGen-calibrated thresholds, one consequence annotation, the filter built into the
  design, shared bootstrap resamples, paired intervals for every gain, dot-and-interval charts.
- The AVI training-chromosome check, now on the primary contrast.

Out:

- The inherited cancer panel, both the ranking test and the PP3 test that issue #2 would keep (Tim's decision). The
  first run's results stay as they are.
- Evo 2. It can be added later on the same indel cases without changing any other number, as an exploratory
  analysis.
- Recessive diagnoses as two different variants (compound heterozygotes); scores rank single variants.
- VEP or any new annotation tool.

## Smallest gain of interest

A gain under 5 points at top 1 would not change what the talk tells a lab to do. In the first run, the locked grpmax
interval for AVI added to the lookups on 500 cases was -3.6 to +3.8 (about plus or minus 3.7 points). At that
precision, a true gain of 5 points gives an interval above zero about three times in four, and a true gain of zero
gives an interval inside -5 to +5 only about half the time. This run's lists will be longer, so intervals may be
wider. A result is reported as "a gain that matters", "no gain that matters" or "inconclusive" by the rule fixed in
the protocol (section 8), not by whether its interval excludes zero.

## What the talk uses

Slides `rd-design`, `rd-results`, `rd-combine`, `rd-indels`, `rd-limits`, `matrix` and `takeaways` will cite this
run's results once it is done. `evee/`, `gpu-drift/`, `evo2-pipeline/` and `brca1-case-ranking/` are not affected.

## Cost

CPU only. Network reads: gnomAD's HGDP + 1KG genotypes at the case positions (about an hour with parallel reads),
CADD by remote tabix, and about 60 MB of EVEE columns. A later Evo 2 run on the indels would be about 2,100 variants,
roughly 30 minutes and USD 3 on one H100 (the first run scored 2,814 indel and cancer variants in 40 minutes for
about USD 3).

## Risks

- Without the family's alleles, more background variants pass the filter, so lists will be longer and top 1 lower
  than in the first run's grpmax analysis (a median of 51 and 34.5 variants, lookups 74 to 80% top 1). How much is
  not known until the frequencies exist. If lists stay short and the lookups stay near 80%, there is little room for
  any gain.
- Homozygous diagnostics in biallelic genes are realistic for consanguineous families but less common in an
  outbred clinic.
- The ClinVar 2-star P/LP diagnostics were often classified partly with PP3 from the same lookups, and the PP3
  thresholds were calibrated on ClinVar. That favours the lookups and both baselines, so Q2 is optimistic and AVI's
  gain is, if anything, conservative.
- AVI's training labels include the 1000 Genomes people; the training-chromosome check measures how much that
  matters.
- The replication pool, after excluding every earlier diagnostic position, may run short in the non-coding class. The
  protocol reduces that quota and logs it rather than reusing positions.
