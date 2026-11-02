# EVEE charts

Two charts on the talk's EVEE slides, redrawn from the figure data that Goodfire published with the EVEE manuscript
([goodfire-ai/evee-manuscript](https://github.com/goodfire-ai/evee-manuscript), `artifacts/*.feather`). No new
analysis: the values are EVEE's, drawn in the talk's colours.

Paper: [paper/paper.pdf](paper/paper.pdf). Visual companion: https://rewire-bio.github.io/amp-2026-studies/evee/

| Output | Slide | From | Shows |
|---|---|---|---|
| `evee_confounding[_gel].svg` | `evee-confounding` | `dataset_characterization.feather` | Share of ClinVar SNVs labelled pathogenic, by consequence, in EVEE's full set of 834,000 variants. Consequence alone predicts much of the label |
| `evee_results[_gel].svg` | `evee-results` | `deconf_benchmark.feather` | AUROC within each consequence on EVEE's deconfounded ClinVar set, for the EVEE probe, Evo 2 zero-shot, AlphaMissense, REVEL and CADD v1.7. Points with fewer than 100 variants are left out |

```sh
# from the repository root, with evee-manuscript cloned next to it
uv run --with matplotlib --with pandas --with pyarrow python evee/evee_charts.py --evee ../evee-manuscript/artifacts
uv run --with matplotlib --with pandas --with pyarrow python evee/evee_charts.py --evee ../evee-manuscript/artifacts --brand gel
```

The SVGs are written to `evee/out/` (not committed). Copy them into the two slide folders in `amp-2026-slides`. Each
chart is drawn for a 1234 × 460 px slot.
