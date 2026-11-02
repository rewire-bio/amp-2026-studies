# Rare-disease case tests, second run

A second run of the rare-disease case tests in [`../rare-disease/`](../rare-disease/), on the intellectual disability
panel only, with frequencies a real patient's pipeline would see. It fixes the problems found in review of the first
run: the 1000 Genomes background people are in gnomAD
([issue #1](https://github.com/rewire-bio/amp-2026-studies/issues/1)), the realistic grpmax filter was only a side
check, the combined model is a weak baseline within each variant type, two consequence annotations were mixed, and
not every gain had a paired interval ([issue #2](https://github.com/rewire-bio/amp-2026-studies/issues/2)).

Status: **protocol written, not yet run.** No results yet. The first run's results are unchanged and still cited at
tag `amp-2026-v4`.

| Document | What it holds |
|---|---|
| [`PROPOSAL.md`](PROPOSAL.md) | Why the rerun, the questions, scope, smallest gain of interest, cost and risks |
| [`PROTOCOL.md`](PROTOCOL.md) | The design, fixed before the 500 replication cases are built, and the deviation log |
| [`PLAN.md`](PLAN.md) | The work, as phases and tasks, with open items and "done when" |

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

## Layout

Folders are created when first used and follow [`../rare-disease/`](../rare-disease/):

```
rare-disease-v2/
  PROPOSAL.md, PROTOCOL.md, PLAN.md
  requirements.txt  pinned package versions
  *.py, *.sh        scripts, copied from rare-disease/ at amp-2026-v4 and changed here, or new
  data/             cases, per-case scores, frequencies (described in data/README.md)
  results/          result files, case checksums, feasibility notes (described in results/README.md)
  figures/          chart script; charts in figures/out/ (not committed)
  logs/             run logs
```

Large inputs (1000 Genomes genotypes for the first run's people, ClinVar, MANE, PanelApp, the AVI table,
AlphaMissense) are read from `../rare-disease/data/`, not copied. REVEL v1.3 is read from
`../../evo2-cloud-operations/data/raw/current-predictors/`, and hg38 from `../../evo2-size-by-use-case/data/raw/hg38.fa.gz`
(decompressed to `data/spliceai/hg38.fa` for SpliceAI). CADD and REVEL scores are not committed, as in the first run
(licence). The EVEE columns used (about 60 MB of the 35 GB Zenodo 19701997 table, CC BY 4.0) are kept in
`data/evee/`, not committed; the per-case EVEE scores are.
