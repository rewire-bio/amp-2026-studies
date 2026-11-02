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
| Frequencies | gnomAD genomes, including the 1000 Genomes samples | gnomAD genomes without the 1000 Genomes samples |
| Frequency filter | None in the main results; AF and grpmax filters as side checks | grpmax below 0.1% built into the cases; grpmax below 1% as a check |
| Consequence | ClinVar MC for groups, CADD for model classes | CADD for SNVs, `annotate_indels.py` for indels, everywhere |
| Baselines | Lookups combined | Rule-based (primary) and lookups combined |
| Added scores | AVI, GPN-MSA, GPN-Star V and M, Evo 2 | AVI (primary); GPN-MSA, GPN-Star V and M (secondary); EVEE in a separate ClinVar-covered analysis; Evo 2 later, if needed |
| Intervals | One random stream, so adding a score moved other intervals | Resamples shared per stratum; a paired interval for every gain |
| Charts | Bars | Dot and interval |

## Layout

Folders are created when first used and follow [`../rare-disease/`](../rare-disease/):

```
rare-disease-v2/
  PROPOSAL.md, PROTOCOL.md, PLAN.md
  *.py, *.sh        scripts, each copied from rare-disease/ at amp-2026-v4 and changed here
  data/             cases, per-case scores, frequencies (described in data/README.md)
  results/          result files, case checksums, feasibility notes (described in results/README.md)
  figures/          chart script; charts in figures/out/ (not committed)
  logs/             run logs
```

Large inputs (1000 Genomes genotypes, ClinVar, MANE, the AVI table, AlphaMissense) are read from
`../rare-disease/data/`, not copied. CADD and REVEL scores are not committed, as in the first run (licence). The EVEE
table (Zenodo 19701997, CC BY 4.0, about 35 GB) stays local in `data/evee/`.
