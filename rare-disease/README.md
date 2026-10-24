# Rare-disease case benchmark

Where does each variant score rank the diagnostic variant in a realistic case? Each case is one 1000 Genomes
high-coverage individual's rare variants (allele frequency below 1%) in a PanelApp panel's green genes, plus one
ClinVar pathogenic or likely pathogenic variant with at least two review stars. For each score we count the
background variants that score equal to or higher than the diagnostic one. Method details are in each script's
docstring; `case_benchmark.py` has the full design.

| Study | Cases | Command | Result |
|---|---|---|---|
| SNVs, panel 285 (intellectual disability) | 200 | `case_benchmark.py` | `results/results.json` |
| Indels up to 50 bp, same panel | see result | `case_benchmark.py --study indel` | `results/results_indel.json` |
| Cancer predisposition, panel 245 | see result | `case_benchmark.py --study cancer` | `results/results_cancer.json` |
| Replication of the SNV study on new cases and backgrounds | 500 | `case_benchmark.py --study replication` | `results/results_replication.json`, [`REPLICATION.md`](REPLICATION.md) |

Scores: CADD v1.7, SpliceAI, REVEL, AlphaMissense (lookups), GPN-MSA, GPN-Star vertebrate and mammal (tables),
and Evo 2 7B scored by us (`evo2/`).

## Rerunning

From this folder, in order (downloads are large; see each script):

```sh
mkdir -p data && cd data && bash ../prepare_data.sh       # panel, MANE, ClinVar, 1000 Genomes
python3 ../build_cases.py [snv|indel|cancer|replication]
python3 ../lookup_scores.py [--study indel]
python3 ../evo2/make_windows.py && bash ../evo2/launch.sh   # optional: Evo 2 on one H100, about USD 20
cd .. && uv run --with pandas --with numpy --with scikit-learn python case_benchmark.py [--study ...]
uv run --with matplotlib python figures/rd_charts.py [--study ...]   # chart SVGs into figures/out/
```

`evo2/run/` and `evo2/run3/` hold the scores and receipts of the two H100 runs (image digest, GPU, timings, cost,
code hash). The raw per-window scores and SSH material are not committed.

## Reproducibility check (24 October 2026)

Rerunning `case_benchmark.py` on the local inputs reproduces `results_indel.json`, `results_cancer.json` and
`results_replication.json` byte for byte. For `results.json`, every row reproduces except "Lookups + Evo 2",
which comes out at top 1 63.0% (52.5 to 73.0) against the committed 62.5% (52.0 to 72.5), and +3.5 points
added to the lookups against the committed +3.0. The Evo 2 score files and the code are unchanged since the
committed run, and scikit-learn 1.7.2 to 1.9.1 give the same rerun value, so the cause is not yet known. The
committed file, which the slides cite, is kept.
