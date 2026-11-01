# GPU drift charts

Charts for the talk's gpu-drift slides, drawn from the recorded BRCA1 runs in
[rewire-bio/evo2-cloud-operations](https://github.com/rewire-bio/evo2-cloud-operations) (clone it next to this
repository). The same 3,893 BRCA1 SNVs were scored with Evo 2 7B on an H100 and on an L40S. With the cut-off fitted on
the H100 scores (the Youden threshold against the Findlay assay labels), 72 calls differ between the two GPUs.

| File | What it shows |
|---|---|
| `flip_charts.py` | Writes the three charts into `out/` (`--brand gel` for the GEL palette) |
| `out/gpu_flip_step1*.png` | One variant's H100 score, just above the cut-off: called loss of function |
| `out/gpu_flip_step2*.png` | The same variant on the L40S, just below it: called functional |
| `out/gpu_flip_scatter*.png` | Every variant, H100 score against L40S score, with the 72 flipped calls marked |

```
uv run --with matplotlib --with pandas python gpu-drift/flip_charts.py [--brand gel]
```
