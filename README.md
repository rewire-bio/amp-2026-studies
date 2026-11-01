# AMP 2026 studies

Code and results behind Tim Richardson's AMP 2026 talk, "Genomic foundation models in the diagnostic lab"
(14 November 2026). The talk itself (slides, narration, companion) is in
[rewire-bio/amp-2026-gfm-talk](https://github.com/rewire-bio/amp-2026-gfm-talk).

All inputs are public. Downloaded data, licensed score tables and run secrets are not committed (see `.gitignore`).

| Folder | Question | Used on slides |
|---|---|---|
| [`rare-disease/`](rare-disease/) | In a case built from a real 1000 Genomes background plus one ClinVar P/LP variant, where does each score rank the diagnostic variant? SNVs, indels, cancer predisposition, and a 500-case replication | `rd-design`, `rd-results`, `rd-combine`, `rd-indels`, `rd-limits`, `cancer-germline`, `cancer-pp3`, `matrix` |
| [`rare-disease/evo2/`](rare-disease/evo2/) | Evo 2 7B zero-shot scores for every case variant, on one Runpod H100, with receipts | `rd-combine`, `scale-cost`, `worked-score` |
| [`brca1-case-ranking/`](brca1-case-ranking/) | On BRCA1 saturation genome editing labels, does Evo 2 put the loss-of-function variant near the top of a simulated case? | backup `bench-*` slides |
| [`evee/`](evee/) | EVEE's published figure data redrawn as two slide charts: pathogenic share by consequence in ClinVar, and AUROC within each consequence | `evee-confounding`, `evee-results` |
| [`evo2-pipeline/`](evo2-pipeline/) | A Nextflow pipeline that adds an Evo 2 score to a case's variants next to CADD and SpliceAI, and ranks the case | `pipeline-overview`, backup `code-*` slides |

Related public repos: [evo2-cloud-operations](https://github.com/rewire-bio/evo2-cloud-operations) (the pinned
Evo 2 scoring image and the BRCA1 study data) and
[goodfire-ai/evee-manuscript](https://github.com/goodfire-ai/evee-manuscript) (EVEE figure data).

The scripts expect sibling checkouts of `evo2-cloud-operations` and `evo2-size-by-use-case` next to this repo for
the REVEL table and the hg38 FASTA.

## Provenance

Moved from the private `rewire-bio/amp-2026-decks` at commit `9a7f2e1`: `rare-disease/` from `v9/rare_disease`,
`v9/evo2` and `v9/figures/rd_charts.py`; `brca1-case-ranking/` from `v4/benchmark` and `v9/benchmark`;
`evo2-pipeline/` from `v3/pipeline`. Paths were updated; the code is otherwise unchanged. The talk's companion
document cites files here at tag `amp-2026-v1`.
