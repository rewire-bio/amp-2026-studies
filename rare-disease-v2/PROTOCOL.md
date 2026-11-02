# Protocol: rare-disease case tests, second run

Status: draft, 2 November 2026. Not yet frozen. Items marked **To confirm** are decided by Tim before the freeze
(section 12). After the freeze, any change goes in the deviation log (section 13) with its date and reason.

Why this run exists and what it asks: [`PROPOSAL.md`](PROPOSAL.md). Work plan: [`PLAN.md`](PLAN.md).

## 1. Frequencies without 1000 Genomes

Every frequency used in this study is a gnomAD genome frequency (v3.1.2, the same 76,215 genomes as v4.1) with the
1000 Genomes samples removed, overall and for each grpmax group (afr, amr, eas, nfe, sas; gnomAD leaves ami, asj,
fin, mid and remaining out of genome grpmax).

The source moves from MyVariant.info, which gives AF only, to the gnomAD release itself (sites VCFs by remote tabix,
or the gnomAD API), because subtraction needs allele counts (AC) and allele numbers (AN).

- **Method A.** If the release carries counts for the 1000 Genomes subset per group, then for each group
  AF = (AC minus AC_1kg) / (AN minus AN_1kg), and grpmax is the highest of the five group AFs.
- **Method B (fallback).** Count alleles in the 1000 Genomes high-coverage VCF (3,202 samples) for the samples that
  are in the gnomAD release, each assigned to the gnomAD group that gnomAD's HGDP + 1KG sample metadata gives it,
  and subtract those counts.

Which method is used is recorded in section 13 before the freeze. A variant with no remaining allele after
subtraction has AF 0 and is kept, as a pipeline would keep a variant it cannot find.

## 2. Frequency filter, built into the design

- **Main filter (dominant mode).** Keep a variant only if its grpmax AF is below 0.1%.
- **Recessive-mode check.** The same cases with grpmax AF below 1%, reported once as a secondary result.
- Diagnostic variants are drawn only from ClinVar candidates that pass the main filter, so every planned case
  exists after filtering and case counts are fixed before any score is looked at.
- Background variants are filtered when the cases are built, not at analysis time.

## 3. One consequence annotation

The same annotation chooses the diagnostic variants, sets the model's classes and defines the reported groups.

- **SNVs:** CADD v1.7's most severe consequence, using the rule `lookup_scores.py` already applies. CADD is looked up
  for the whole ClinVar candidate pool before cases are built.
- **Indels:** `annotate_indels.py` (MANE Select).

Classes:

| Class | SNV (CADD) | Indel (`annotate_indels.py`) |
|---|---|---|
| nonsense / lof | `STOP_GAINED`, `STOP_LOST` | |
| frameshift | | frameshift |
| inframe | | inframe |
| missense | `NON_SYNONYMOUS` | |
| canonical splice | `CANONICAL_SPLICE` | canonical splice |
| splice region | `SPLICE_SITE` | |
| synonymous | `SYNONYMOUS` | |
| UTR | `5PRIME_UTR`, `3PRIME_UTR` | UTR |
| intronic | `INTRONIC` | intronic |
| other | anything else | anything else |

Canonical splice and splice region are kept apart; the first run merged them in the model.

## 4. Cases

| Study | Cases | People | Diagnostic quotas |
|---|---|---|---|
| SNV | 200, 4 per person | the first run's 50 (`../rare-disease/data/samples50.txt`) | missense 60, canonical splice 40, nonsense 30, other 70 (synonymous, UTR, intronic, splice region) |
| Indel (up to 50 bp) | 120 | the same 50 | frameshift 40, inframe 30, canonical splice 30, UTR or intronic 20 |
| Replication | 500, 10 per person | the first run's 50 replication people (`../rare-disease/data/replication/samples50.txt`) | missense 150, canonical splice 100, nonsense 75, other 175; no position used as an SNV diagnostic |

- Inputs fixed to the first run's: ClinVar fileDate 2026-10-11, P/LP with 2+ review stars; PanelApp 285 v11.45 green
  genes; MANE v1.5 Select exons plus 20 bp; 1000 Genomes high-coverage phased panel (20220422).
- Background: the person's variants with 1000 Genomes AF below 1% in the regions, then the main filter. Background
  variants that are ClinVar P/LP are removed, so each case has one right answer.
- Seeds: SNV 20261102, indel 20261103, replication 20261109.
- The replication cases are built only after the freeze (section 12).

## 5. Scores

| Score | Role | Covers |
|---|---|---|
| CADD v1.7 PHRED | lookup | all SNVs; indels seen in gnomAD |
| SpliceAI max delta | lookup | SNVs from CADD's annotation; indels run ourselves (`run_spliceai.sh`) |
| REVEL v1.3 | lookup | missense |
| AlphaMissense | lookup | missense |
| AlphaGenome AVI | primary added score | SNVs |
| GPN-MSA | secondary added score | SNVs |
| GPN-Star V and M | secondary added scores | SNVs |
| EVEE | secondary, own analysis (section 9) | ClinVar variants only |

- A missing score ranks last, as a failed lookup would.
- CADD and precomputed SpliceAI stay unavailable for indels in the models, as in the first run: their indel
  coverage comes from gnomAD and would give the answer away.
- AVI is the one primary added score. GPN-MSA, GPN-Star V, GPN-Star M and EVEE are secondary. Results state that
  there are four secondary comparisons, so a single interval above zero among them is not read as a finding.

## 6. Baselines

- **Lookups combined.** Logistic regression on the four lookups, consequence class and score-by-class terms, as in
  the first run.
- **Rule-based (primary comparator).** Each variant is scored by the tool for its class: SpliceAI for canonical
  splice and splice region, the higher of REVEL and AlphaMissense for missense, CADD for the rest. The scores are put
  on one scale by their calibrated PP3 evidence level (none, supporting, moderate, strong), using thresholds from
  Pejaver et al. 2022 (REVEL, CADD), Bergquist et al. 2025 (AlphaMissense, REVEL) and Walker et al. 2023 (SpliceAI).
  Ties within a level are broken by CADD. The thresholds are copied from the papers into this section, with table
  numbers, before the freeze. **To confirm:** the shared scale and the tie-break.

Every gain is reported against both baselines. The rule-based baseline is primary because it matches the best
single tool for each type.

## 7. Models and test sets

- Cross-validated: 5-fold, grouped by person, on the 200 SNV cases and on the 120 indel cases.
- Locked: every model fitted once on all 200 SNV cases and applied unchanged to the 500.
- Combinations: lookups combined; lookups plus one added score at a time. Adding a score never changes another
  combination.

## 8. Outcomes

- **Primary.** Top-1 gain of lookups plus AVI over the rule-based baseline, locked, on the 500 cases, with the main
  filter. The same gain over lookups combined is reported alongside.
- **Secondary.** Top 5 and top 10; each class with at least 20 cases; the 200 cases (cross-validated); the indel
  cases; the recessive-mode check; the secondary scores; the AVI held-out chromosome check (`avi_heldout.py` logic:
  AVI's training chromosomes are 1, 4, 7, 8, 10, 13 and 15).
- Smallest gain of interest: 5 points at top 1 ([`PROPOSAL.md`](PROPOSAL.md)).

## 9. EVEE: separate analysis

EVEE's released table (Zenodo 19701997, CC BY 4.0) has one row per ClinVar variant (ClinVar 2026-03-21). Every
diagnostic variant is in ClinVar and many background variants are not, so in the main design EVEE would rank the
diagnostic first by coverage alone. Its probe was also trained on ClinVar P/LP against B/LB labels. So EVEE is never
in the main results.

- **Feasibility, before the freeze.** Measure how many background variants per case have an EVEE row. Check the
  manuscript and the Zenodo record for whether the released scores are cross-validated (held-out) predictions.
- **Rule.** If the median case keeps fewer than 20 background variants with an EVEE row, EVEE is reported as not
  testable here, with the coverage figures. **To confirm:** the threshold of 20.
- **Otherwise.** A "ClinVar-covered" analysis: every case's list is restricted to variants with an EVEE row, and every
  score is ranked on those same lists. Where numbers allow, diagnostic variants are limited to those not P/LP in the
  2026-03-21 ClinVar release, so they cannot be EVEE training labels. Results are labelled as not comparable with
  the main analysis.

## 10. Adding Evo 2 later

The design must let Evo 2 be added without changing any other number.

- The case files are frozen: SHA-256 checksums of every `cases.tsv` and `case_variants.tsv` are written to
  `results/case_checksums.txt` and tagged.
- `evo2/make_windows.py` takes a required `--fasta`; the README says where hg38 comes from.
- Intervals use shared resamples (section 11), so adding a score moves no other interval.
- Rule for running Evo 2 on the frozen indel cases: run it if lookups combined put the diagnostic indel first in
  fewer than 80% of indel cases. **To confirm:** the 80% threshold.
- Check after adding it: every row of every result file that does not involve Evo 2 is byte-for-byte the same.

## 11. Intervals

- 1,000 bootstrap resamples of people (each with all their cases).
- Drawn once per stratum, from a seed derived from the stratum name, and reused by every score and every paired
  difference in that stratum.
- Every gain has a paired 95% interval.
- Charts are dot-and-interval, not bars.

## 12. Freeze

Before the 500 replication cases are built, this protocol, the analysis code and the 200 and 120 case files are
committed and tagged `rd-v2-protocol`. The open items marked **To confirm** must be closed by then:

- Method A or B for frequencies (section 1).
- The rule-based baseline's scale, thresholds and tie-break (section 6).
- EVEE go or not testable (section 9).
- The Evo 2 threshold (section 10).

## 13. Deviation log

| Date | Section | Change | Reason |
|---|---|---|---|
| | | | |
