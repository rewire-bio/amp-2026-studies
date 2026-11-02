# Rare-disease case tests, second run

A second run of the rare-disease case tests in [`../rare-disease/`](../rare-disease/), on the intellectual disability
panel only, with frequencies a real patient's pipeline would see. It fixes the problems found in review of the first
run: the 1000 Genomes background people are in gnomAD
([issue #1](https://github.com/rewire-bio/amp-2026-studies/issues/1)), the realistic grpmax filter was only a side
check, the combined model is a weak baseline within each variant type, two consequence annotations were mixed, and
not every gain had a paired interval ([issue #2](https://github.com/rewire-bio/amp-2026-studies/issues/2)).

Status: **done** (2 November 2026). The results are the files at the current commit. No release tag yet; the talk
still cites the first run at tag `amp-2026-v4`. Evo 2 on the indel cases (PLAN.md phase 6) has not been run.

| Document | What it holds |
|---|---|
| [`PROPOSAL.md`](PROPOSAL.md) | Why the rerun, the questions, scope, smallest gain of interest, cost and risks |
| [`PROTOCOL.md`](PROTOCOL.md) | The design and the deviation log |
| [`PLAN.md`](PLAN.md) | The work, as phases and tasks |
| [`results/analysis_replication.md`](results/analysis_replication.md) | The 500 replication cases, locked: the primary result and every check |
| [`results/analysis_snv_indel.md`](results/analysis_snv_indel.md) | The 200 SNV and 120 indel cases, cross-validated (descriptive) |
| [`results/README.md`](results/README.md) | Every result file |
| [`data/README.md`](data/README.md) | Every data file, and what is not committed |

## Results

All figures are top 1 (the share of cases where the diagnostic variant ranks first) unless stated, in percent, with
95% bootstrap intervals over people. Gains are in percentage points with paired intervals. The 500 replication cases
use models fitted on the 200 and locked. Main filters: grpmax AF below 0.1% with the person and their family left
out of gnomAD, and the background inheritance filter.

### Q1 (primary): does adding AVI to the lookups raise top 1 on the 500?

No gain that matters. Adding AlphaGenome AVI to the lookups combined raises top 1 from 90.6% to 92.4%, a gain of
+1.8 points (-0.4 to +4.0). The interval lies inside -5 to +5, the reading rule fixed in advance for "no gain that
matters". The lookups already rank the diagnostic first in 9 cases in 10, and lists are short (median 19 variants).

| 500 cases, locked | Without AVI | With AVI | Gain |
|---|---|---|---|
| Lookups combined | 90.6 (87.8 to 93.4) | 92.4 (89.8 to 94.8) | +1.8 (-0.4 to +4.0) |
| Rule-based | 87.2 (82.8 to 91.0) | 89.8 (86.6 to 92.8) | +2.6 (+0.4 to +5.6) |

On the 200 cases, cross-validated and descriptive, the primary gain is +4.0 (-2.0 to +9.5).

### Q2: how often do the lookups alone rank the diagnostic first?

In 9 cases in 10. On the 500, the lookups combined by the locked model reach 90.6% (87.8 to 93.4) and the rule-based
ranking 87.2% (82.8 to 91.0); the difference is +3.4 (0.0 to +7.4). At top 5 they reach 97.0% (95.4 to 98.4) and
97.8% (96.6 to 98.8).

| 500 cases | Cases | Lookups combined (locked) | Rule-based | Best lookup |
|---|---|---|---|---|
| Canonical splice | 100 | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | SpliceAI 96.0 (91.1 to 100.0) |
| Nonsense | 75 | 100.0 (100.0 to 100.0) | 100.0 (100.0 to 100.0) | CADD 98.7 (95.5 to 100.0) |
| Missense | 150 | 87.3 (82.1 to 92.5) | 90.0 (84.0 to 95.4) | REVEL 90.7 (84.7 to 95.7) |
| Non-coding | 175 | 84.0 (77.4 to 89.9) | 72.0 (61.7 to 81.1) | SpliceAI 86.3 (81.3 to 91.3) |

Lookups minus rule is +12.0 (+4.7 to +20.9) on non-coding diagnostics, where the rule has just SpliceAI's supporting
level to go on. On missense the rule is ahead on the estimate (-2.7, -7.3 to +2.8); on canonical splice and nonsense
the two are equal. Within each type the combined model is now close to the best lookup, unlike the first
run (issue #2, problem 3). On the 200 SNV cases (cross-validated) the lookups reach 84.0% (75.5 to 91.5) and the rule 81.0%
(73.5 to 88.0). On the 120 indel cases, ranked among a median of 4 indels, they reach 84.2% (75.4 to 92.4) and
76.7% (65.8 to 86.8).

### Q3: does AVI help the rule-based ranking too?

A little. Rule + AVI minus Rule-based is +2.6 (+0.4 to +5.6) on the 500: above zero but under the 5-point smallest
gain of interest, so it does not change the answer to Q1. The gain is in non-coding diagnostics (+8.6, +3.2 to
+15.1), where AVI separates variants the rule scores alike. On the 200 the key secondary is +4.5 (+0.5 to +8.0), but
it sits in the cases whose diagnostic is absent from gnomAD and is flagged as possibly due to gnomAD presence in the
all1kg check; on the 500 it is not.

### Q4: do GPN-MSA or GPN-Star add anything?

No. Added to the lookups on the 500, each with a Bonferroni-adjusted 98.33% interval: GPN-MSA -0.2 (-1.6 to +1.4),
GPN-Star-V +0.4 (-1.8 to +2.4), GPN-Star-M -0.2 (-1.6 to +1.2). None excludes zero.

### EVEE (descriptive)

On the 500 with the 1% filter, lists restricted to variants with an EVEE row (492 cases), EVEE alone ranks the
diagnostic first in 87.0% (82.5 to 90.9) of cases, the same as the locked lookups, 87.0% (82.2 to 91.0). This is
agreement with ClinVar labels on in-sample variants, not diagnostic performance: EVEE covers ClinVar variants only,
and almost every diagnostic was already P/LP in its training labels.

### Checks

| Check (500, locked) | Lists, median | Lookups combined | AVI gain to the lookups | AVI gain to the rule |
|---|---|---|---|---|
| Main | 19 | 90.6 (87.8 to 93.4) | +1.8 (-0.4 to +4.0) | +2.6 (+0.4 to +5.6) |
| Background grpmax below 1% | 69.5 | 82.0 (76.4 to 87.2) | -1.2 (-6.6 to +3.2) | +5.2 (+1.4 to +9.6) |
| No background inheritance filter | 34.5 | 79.6 (73.2 to 85.4) | +4.0 (0.0 to +8.0) | +2.2 (-1.2 to +5.4) |
| All of 1000 Genomes left out | 21 | 87.8 (83.2 to 91.8) | +0.4 (-2.8 to +3.0) | +2.6 (+0.4 to +5.4) |

- **AVI training chromosomes.** The primary gain is no larger on the chromosomes AVI was trained on: training minus
  other -0.5 (-3.7 to +3.1).
- **Diagnostic in gnomAD or not.** 114 of the 500 diagnostics are in gnomAD (AC above 0) and 386 are not. The primary
  gain is +2.6 (+0.5 to +4.8) in absent cases and -0.9 (-6.6 to +5.0) in found ones; absent minus found is +3.5 (-2.2
  to +8.8), so it is not flagged. Neither is the key secondary.
- **Rule ties.** No diagnostic is tied with a background variant under the main filters.
- **Invariance.** Dropping an added score changes no row that does not involve it (`tests/test_invariance.py`).

Details: [`results/analysis_replication.md`](results/analysis_replication.md) and
[`results/analysis_snv_indel.md`](results/analysis_snv_indel.md).

## What this means for the talk

The first run reported that AVI's gain held on the 500 new cases: +6.2 points locked (0.0 to +13.0), and +7.4 (+1.2
to +14.2) after a gnomAD AF filter of 0.1%. Under its stricter grpmax filter the same gain was +0.2 (-3.6 to +3.8).
This run, with frequencies a patient's pipeline would see, a grpmax filter and an inheritance filter, finds +1.8
(-0.4 to +4.0): no gain that matters. It agrees with the first run's grpmax result, not its headline.

So the talk should not say that AVI's gain holds on new cases. It can say that with a realistic pipeline the
established lookups rank the diagnostic first in about 9 cases in 10, that a simple rule-based ranking does nearly as
well, and that AVI helps a little with non-coding diagnostics only. GPN-MSA and GPN-Star add nothing.

## Limits

- **Short lists.** After the filters the median list is 19 SNVs on the 500 and 4 indels on the 120. A real exome or
  genome has more candidates; these lists hold only rare variants in the panel's genes.
- **Lookups near the ceiling.** At 90.6% top 1 (87.8 to 93.4) there is little room for any added score to gain.
- **Homozygous diagnostics.** Diagnostics in biallelic genes are homozygous (293 of the 500,
  [`results/build_replication.md`](results/build_replication.md)). That fits a consanguineous family; in an outbred
  clinic most recessive diagnoses are two different variants, which this design does not model.
- **ClinVar diagnostics were partly classified with the same lookups.** Many 2-star P/LP variants used PP3 evidence
  from these tools, and the PP3 thresholds were calibrated on ClinVar. That favours the lookups and both baselines.
- **1000 Genomes people are not patients.** Their backgrounds hold no disease variant of their own, and the
  diagnostic is inserted, not observed.
- **gnomAD presence is only partly fixed.** Diagnostics keep their full gnomAD frequency by design, and most are
  absent from gnomAD (AC 0): 148 of the 200 (74.0%) and 386 of the 500 (77.2%). Background variants rarely are: on the
  first run's 200 SNV cases, 2% once the family is left out
  ([`results/gnomad_overlap.md`](results/gnomad_overlap.md)). The found and absent split (above) shows no gain driven
  by it on the 500.
- **Cross-validated models move a few points when folds change.** People are assigned to folds by row count, so a
  small change to the lists moves people between folds. The 200 and 120 results are descriptive for this reason.
- **Indel lists are short.** With a median of 4 indels per case, the indel results say little about how scores
  separate many candidates.

## What changes from the first run

| | First run (`rare-disease/`) | This run |
|---|---|---|
| Panels | 285 intellectual disability; 245 inherited cancer | 285 only |
| Frequencies | gnomAD genomes, including the person and all 1000 Genomes | gnomAD genomes without the person and their family; without all of 1000 Genomes as a check |
| Frequency filter | None in the main results; AF and grpmax filters as side checks | grpmax below 0.1% in the main analysis; below 1% as a check |
| Inheritance | Ignored: many diagnostics were carrier states | Every diagnostic fits its gene's inheritance; background filtered by inheritance |
| Consequence | ClinVar MC for groups, CADD for model classes | CADD (MANE transcripts) for SNVs, `annotate_indels.py` for indels, everywhere |
| Baselines | Lookups combined | Lookups combined, and a rule-based ranking from ClinGen-calibrated thresholds |
| Primary contrast | Lookups plus AVI against lookups combined | The same, with a reading rule fixed in advance and AVI also tested on the rule-based ranking |
| Replication people | 50 new people | 50 people used in neither of the first run's sets |
| Added scores | AVI, GPN-MSA, GPN-Star V and M, Evo 2 | AVI (primary); GPN-MSA, GPN-Star V and M (secondary, adjusted intervals); EVEE descriptive only; Evo 2 later, exploratory |
| Intervals | One random stream, so adding a score moved other intervals | Resamples shared per stratum; a paired interval for every gain |
| Charts | Bars | Dot and interval |

## How to rerun

Run from `rare-disease-v2/`. `RUN` stands for `uv run -q --python 3.12 --with-requirements requirements.txt python`.
Each script's docstring has the details. The first run's inputs must be in `../rare-disease/data/` (see its README),
REVEL in `../../evo2-cloud-operations/data/raw/current-predictors/`, and hg38 in
`../../evo2-size-by-use-case/data/raw/hg38.fa.gz`. The shell scripts need `bcftools`, `tabix`, `samtools` and `curl`.

```
RUN people.py                              # data/people.tsv, data/rep/samples50.txt
bash fetch_rep_kgp.sh                      # data/rep/kgp/chrN.bcf
RUN backgrounds.py main rep                # data/<set>/person_variants.tsv
RUN gnomad_af.py full                      # data/gnomad_full.tsv
RUN gnomad_af.py fetch                     # data/gnomad_subset/dense/ (about an hour)
RUN gnomad_af.py counts                    # data/gnomad_subset/counts_chrN.tsv
RUN gnomad_af.py af                        # data/<set>/person_variant_af.tsv
RUN gnomad_af.py checks                    # results/gnomad_checks.md
RUN gnomad_af.py overlap                   # results/gnomad_overlap.md
RUN moi.py                                 # data/gene_moi.tsv
RUN pool.py                                # data/pool/, results/pool_counts.md
RUN build_cases.py --study snv             # then --study indel and --study replication
RUN annotate_indels.py --study snv         # then --study indel
gunzip -c ../../evo2-size-by-use-case/data/raw/hg38.fa.gz > data/spliceai/hg38.fa && samtools faidx data/spliceai/hg38.fa
bash run_spliceai.sh                       # data/spliceai/indels.out.vcf
RUN lookup_scores.py --study snv           # then --study indel and --study replication
RUN case_benchmark.py --study snv --analysis main
                                           # every study with main, bg1pct, nomoi, all1kg;
                                           # snv and replication also with heldout and evee
RUN figures/rd_charts.py                   # figures/out/ (not committed)
RUN tests/test_invariance.py replication   # and tests/test_rule_ties.py
```

`gnomad_af.py all` runs its six steps in order. CADD and REVEL scores are not committed (licence), so
`lookup_scores.py` fetches them again; CADD needs no download (remote tabix).

## Layout

```
rare-disease-v2/
  PROPOSAL.md, PROTOCOL.md, PLAN.md
  requirements.txt  pinned package versions
  *.py, *.sh        scripts, copied from rare-disease/ at amp-2026-v4 and changed here, or new
  data/             cases, per-case scores, frequencies (described in data/README.md)
  results/          result files and notes (described in results/README.md)
  figures/          chart script; charts in figures/out/ (not committed)
  tests/            invariance and rule-tie tests
  logs/             run logs
```

Large inputs (1000 Genomes genotypes for the first run's people, ClinVar, MANE, PanelApp, the AVI table,
AlphaMissense) are read from `../rare-disease/data/`, not copied. REVEL v1.3 is read from
`../../evo2-cloud-operations/data/raw/current-predictors/`, and hg38 from `../../evo2-size-by-use-case/data/raw/hg38.fa.gz`
(decompressed to `data/spliceai/hg38.fa` for SpliceAI). What is not committed, and why, is listed in
[`data/README.md`](data/README.md).
