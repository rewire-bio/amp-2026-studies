# Research proposal: rare-disease case tests, second run

Status: proposed, 2 November 2026. Design in [`PROTOCOL.md`](PROTOCOL.md), work plan in [`PLAN.md`](PLAN.md).
Tracked in issues [#1](https://github.com/rewire-bio/amp-2026-studies/issues/1) and
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

A review of the talk found five problems with that design.

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
   matching single tool in every type: canonical splice 82% against SpliceAI 90%, nonsense 60% against CADD 90%,
   missense 43% against REVEL or AlphaMissense 58%. Part of an added score's gain may make up for that.
4. **Two consequence annotations.** Diagnostic variants were grouped by ClinVar's MC field; the model's classes came
   from CADD, whose splice class also takes splice-region variants.
5. **Small groups.** With 30 to 70 cases per type, most per-type differences are noise, and not every gain had a
   paired interval.

## Questions

- **Q1 (primary).** With frequencies a real patient's pipeline would see (gnomAD without the 1000 Genomes samples)
  and a grpmax filter below 0.1%, does adding AVI to the lookups raise top-1 on 500 new cases, against a
  rule-based baseline and against the combined lookups?
- **Q2.** How often do the lookups alone rank the diagnostic variant first, overall and by variant type?
- **Q3.** On indels, how much room is left after the lookups? This decides whether Evo 2 is run on the indel cases.
- **Q4 (secondary).** Do GPN-MSA, GPN-Star (vertebrate and mammal) or EVEE add anything under the same design?

## Scope

In:

- PanelApp 285 (intellectual disability) only.
- 200 SNV cases and 120 indel cases from the first run's 50 people; 500 SNV replication cases from the other 50.
- Scores: CADD v1.7, SpliceAI, REVEL, AlphaMissense (the lookups), AlphaGenome AVI (primary added score), GPN-MSA
  and GPN-Star V and M (secondary), EVEE (secondary, in a separate analysis because it covers ClinVar variants only;
  see the protocol).
- A rule-based baseline, one consequence annotation, the filter built into the case design, paired intervals for
  every gain, dot-and-interval charts.
- The AVI held-out chromosome check.

Out:

- The inherited cancer panel (both the ranking test and the PP3 test). The first run's results stay as they are.
- Evo 2. It can be added later on the frozen indel cases without changing any other number (protocol section 10).
- VEP or any new annotation tool.

## Smallest gain of interest

A gain under 5 points at top 1 would not change what the talk tells a lab to do. On the first run, the locked grpmax
interval for AVI on the 500 cases was about plus or minus 4 points, so 500 cases can tell a 5 point gain from none.
A result whose interval sits inside -5 to +5 will be reported as "no gain that matters", not as "underpowered".

## What the talk uses

Slides `rd-design`, `rd-results`, `rd-combine`, `rd-indels`, `rd-limits`, `matrix` and `takeaways` will cite this
run's results once it is tagged. `evee/`, `gpu-drift/`, `evo2-pipeline/` and `brca1-case-ranking/` are not affected.

## Cost

CPU only, plus downloads: the EVEE table (about 35 GB) and gnomAD sites data for the panel regions. A later Evo 2
run on the indel cases would be a small job on one H100 (the first run scored 13,345 variants in 2 h 50 min for
about USD 11).

## Risks

- gnomAD may not publish 1000 Genomes subset counts that can be subtracted. The protocol has a fallback that counts
  the alleles from the 1000 Genomes genotypes directly.
- After the stricter filter, lists are short (a median of 35 to 51 variants in the first run). Top-1 rates near 80%
  leave little room for any gain.
- The ClinVar pool, after the filter and after excluding the 200's positions, may not fill the 500's quotas for
  some types. The plan checks this before building cases.
- EVEE may cover too few background variants to be tested at all, and its released scores may include ClinVar
  training labels. The protocol says how that is decided.
