# Rare-disease results

Every result of the case benchmark. `case_benchmark.py` writes all files except `avi_heldout.json` and
`compare_amp-2026-v1.txt`. The talk's companion reads values from these files at tag `amp-2026-v3` and checks each
number it quotes against them.

## Files

| File | Study | Filter |
|---|---|---|
| `results.json` | 200 SNV cases | none |
| `results_indel.json` | 120 indel cases | none |
| `results_cancer.json` | 140 inherited cancer cases | none |
| `results_replication.json` | 500 new SNV cases, with locked models | none |
| `results_gnomad0.001.json`, `results_replication_gnomad0.001.json`, `results_indel_gnomad0.001.json`, `results_cancer_gnomad0.001.json` | as above | gnomAD genome AF below 0.1% |
| `results_gnomad0.01.json`, `results_replication_gnomad0.01.json` | SNV, replication | AF below 1% |
| `results_gnomad0.001_grpmax.json`, `results_replication_gnomad0.001_grpmax.json` | SNV, replication | grpmax AF below 0.1% |
| `avi_heldout.json` | SNV, replication | none; AVI's training chromosomes against the rest (`avi_heldout.py`) |
| `compare_amp-2026-v1.txt` | all | what changed when AVI was added (`compare_results.py`) |

## Fields in `results*.json`

```
{
  "cases": 200,                          number of cases
  "median_variants_per_case": 234.5,     background plus diagnostic
  "coverage_all_variants": {score: n},   case-variant rows with a value for each score
  "coverage_diagnostic":   {score: n},   diagnostic variants with a value for each score
  "pp3_missense": {                      ACMG PP3 thresholds, missense variants only
    "AlphaMissense": {"diagnostic": {"n", "supporting", "moderate", "strong"},
                      "background": {...}},       share of variants at or above each threshold
    "REVEL": {...}
  },
  "strata": {
    "all" | <consequence group>: {      all cases, then cases split by diagnostic consequence
      "cases": n,
      "scores": {
        <score or model>: {
          "median_competitors": 4.0,     background variants scoring at least as high as the diagnostic one
          "top_1": 0.28, "top_1_ci95": [low, high],
          "top_5": ..., "top_5_ci95": ...,
          "top_10": ..., "top_10_ci95": ...
        }
      },
      "added_to_lookups": {             paired difference from "Lookups combined", same cases
        "Lookups + <score>": {"top_1", "top_1_ci95", "top_5", ..., "top_10", ...}
      },
      "added_to_locked_lookups": {...}  replication only: the same for the locked models
    }
  },
  "gnomad_filter": {                     filtered runs only
    "field": "AF" | "AF_grpmax", "max_af": 0.001,
    "cases_before": 200, "variants_before": 61064, "variants_after": 23896
  }
}
```

Shares are fractions (0.28 is 28%). Intervals are 95% bootstrap intervals over the 50 background people.

Score and model names:

- **Single scores:** `CADD v1.7`, `SpliceAI`, `REVEL`, `AlphaMissense`, `GPN-MSA`, `GPN-Star-V`, `GPN-Star-M`,
  `AlphaGenome AVI` and `Evo 2 7B`. Evo 2 is not in the replication files.
- **Combined models:** `Lookups combined`, and `Lookups + <score>` for each foundation-model score.
- **Locked models:** in the replication files, the combined models also appear with ` (locked)`. These were fitted
  on the 200 SNV cases and applied unchanged.

Consequence groups (`strata` keys): SNV and replication studies use `missense`, `nonsense`, `canonical_splice` and
`noncoding_or_synonymous`. Indels use `frameshift`, `inframe`, `canonical_splice` and `noncoding`. The cancer study
uses `missense`, `nonsense`, `frameshift`, `canonical_splice` and `noncoding`.

## `avi_heldout.json`

For `snv` and `replication`, two splits: `training_chromosomes` (AVI's training chromosomes 1, 4, 7, 8, 10, 13 and
15) and `other_chromosomes`. Each split has `cases`, `median_variants_per_case`, and `top_1` and `top_5` for every
single score. Cases are assigned by their diagnostic variant's chromosome, and ranked only among background variants
in the same split. Background sizes differ between the splits, so compare scores within a split.

## Reading a number

To find "AVI added to the lookups, 500 new cases, locked, top 1":

```sh
python3 -c "import json; d=json.load(open('results_replication.json')); \
print(d['strata']['all']['added_to_locked_lookups']['Lookups + AlphaGenome AVI (locked)'])"
```
