"""Charts for the gpu-drift paper and companion.

Reads the committed evidence in a sibling evo2-cloud-operations checkout (commit f8ae11c9): the recorded BRCA1 scores,
evidence/classification-flips/results.json and evidence/followup-20261011/followup-results.json.

  uv run --with matplotlib --with pandas python gpu-drift/paper/make_figures.py

Writes PDF into gpu-drift/paper/figures/ and PNG into gpu-drift/figures/. The two flip charts reuse the drawing
functions in ../flip_charts.py unchanged (rewire palette); the two summary charts are drawn here.
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY))
import flip_charts as fc  # noqa: E402

PDF = HERE / "figures"
PNG = STUDY / "figures"
FLIPS = fc.OPS / "classification-flips" / "results.json"
FOLLOWUP = fc.OPS / "followup-20261011" / "followup-results.json"
C = fc.BRANDS["rewire"]
THRESHOLDS = [("youden", "Youden"), ("sensitivity_0.90", "90% sensitivity"), ("sensitivity_0.95", "95% sensitivity")]


def save(fig, name: str) -> None:
    fig.savefig(PDF / f"{name}.pdf", facecolor=C["surface"])
    fig.savefig(PNG / f"{name}.png", dpi=150, facecolor=C["surface"])
    plt.close(fig)


def flips_by_threshold(r: dict) -> None:
    """Share of calls that change, with exact 95% intervals, for each pair and threshold."""
    full = r["full_set"]["recorded: H100 (P2) vs L40S (P5)"]
    sub = r["subset"]
    rows = [
        ("H100 vs L40S, full set (n = 3,893)", full, C["flip"]),
        ("H100 vs L40S, subset, default", sub["default: H100 vs L40S"]["flips"], C["flip"]),
        ("H100 vs L40S, subset, FP8 off", sub["fp8off: H100 vs L40S"]["flips"], C["flip"]),
        ("H100 vs L40S, subset, FP32", sub["fp32: H100 vs L40S"]["flips"], C["dot"]),
        ("L40 vs L40S, subset, default", sub["default: L40 vs L40S"]["flips"], C["flip"]),
        ("L40 vs L40S, subset, FP8 off", sub["fp8off: L40 vs L40S"]["flips"], C["flip"]),
        ("H100 vs H100, two clouds (n = 3,893)", None, C["dot"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 4.8), sharey=True)
    for ax, (key, title) in zip(axes, THRESHOLDS):
        for i, (label, res, col) in enumerate(rows):
            y = len(rows) - 1 - i
            if res is None:  # control: 0 flips at every threshold; exact (Clopper-Pearson) interval for 0 of n
                k, n = r["controls"]["recorded: P2 vs P3 (both H100)"], 3893
                assert k == 0, k
                p, lo, hi = 0.0, 0.0, 1 - 0.025 ** (1 / n)
            else:
                s = res[key]["strata"]["all"]
                k, n, p, (lo, hi) = s["flips"], s["n"], s["proportion"], s["ci95"]
            ax.plot([100 * lo, 100 * hi], [y, y], color=col, lw=2, solid_capstyle="round")
            ax.scatter([100 * p], [y], s=40, color=col, zorder=3, edgecolors=C["surface"], linewidths=1.5)
            ax.text(100 * max(hi, p) + 0.25, y, f"{k}/{n:,}", va="center", fontsize=10, color=C["muted"])
        ax.set_title(title, fontsize=12.5, color=C["dark"])
        ax.set_xlim(0, 10.5)
        ax.set_xlabel("Calls that change (%)", fontsize=11, color=C["dark"])
        ax.tick_params(labelsize=10, colors=C["muted"])
        ax.grid(axis="x", color="#E4E2DE", lw=0.6)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_facecolor(C["surface"])
    axes[0].set_yticks(range(len(rows)))
    axes[0].set_yticklabels([r[0] for r in reversed(rows)], fontsize=10.5, color=C["dark"])
    fig.tight_layout()
    save(fig, "flips_by_threshold")


def differences(f: dict) -> None:
    """Median absolute per-variant score difference for each comparison, log scale."""
    cmp = f["comparisons"]
    rows = [
        ("H100 vs L40S, default", cmp["default: H100 vs L40S"], C["flip"]),
        ("L40 vs L40S, default", cmp["same architecture: second Ada GPU vs L40S, default"], C["flip"]),
        ("H100 vs L40S, FP8 off", cmp["fp8off: H100 vs L40S"], C["flip"]),
        ("L40 vs L40S, FP8 off", cmp["same architecture: second Ada GPU vs L40S, fp8off"], C["flip"]),
        ("H100 vs L40S, FP32", cmp["fp32: H100 vs L40S"], C["dot"]),
        ("H100 only: FP8 off vs default", cmp["H100: fp8off vs default"], C["muted"]),
        ("H100 only: FP32 vs default", cmp["H100: fp32 vs default"], C["muted"]),
        ("L40S only: FP8 off vs default", cmp["L40S: fp8off vs default"], C["muted"]),
        ("L40S only: FP32 vs default", cmp["L40S: fp32 vs default"], C["muted"]),
        ("H100 only: batch 8 vs batch 1", cmp["H100: batch 8 vs batch 1"], C["dot"]),
    ]
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    floor = 1e-7
    for i, (label, c, col) in enumerate(rows):
        y = len(rows) - 1 - i
        med, mx = c["median_abs_difference"], c["max_abs_difference"]
        if med == 0:
            ax.text(floor * 1.3, y, "0 (identical)", va="center", fontsize=10.5, color=C["dark"])
            continue
        ax.plot([med, mx], [y, y], color=col, lw=1.4, alpha=0.5)
        ax.scatter([med], [y], s=42, color=col, zorder=3, edgecolors=C["surface"], linewidths=1.5)
        ax.scatter([mx], [y], s=22, marker="|", color=col, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(floor, 3e-3)
    ax.set_ylim(-0.7, len(rows) - 0.4)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in reversed(rows)], fontsize=10.5, color=C["dark"])
    ax.set_xlabel("Absolute score difference (dot: median; bar end: maximum)",
                  fontsize=11, color=C["dark"])
    ax.tick_params(labelsize=10, colors=C["muted"])
    ax.grid(axis="x", color="#E4E2DE", lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor(C["surface"])
    fig.tight_layout()
    save(fig, "differences")


def main() -> None:
    PDF.mkdir(exist_ok=True)
    PNG.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
    d, cut = fc.scores()
    flips = d[d.flip & (d.h100 >= cut)].assign(gap=lambda x: x.h100 - x.l40s).sort_values("gap", ascending=False)
    v = flips.iloc[0]
    for name, fn in [("gpu_flip_step2", lambda p: fc.number_line(C, p, v, cut, True)),
                     ("gpu_flip_scatter", lambda p: fc.scatter(C, p, d, cut))]:
        fn(PDF / f"{name}.pdf")
        fn(PNG / f"{name}.png")
    flips_by_threshold(json.loads(FLIPS.read_text()))
    differences(json.loads(FOLLOWUP.read_text()))
    print(f"cut-off {cut:.6f}; flips {int(d.flip.sum())} of {len(d)}; example {v.id}: H100 {v.h100:.6f}, "
          f"L40S {v.l40s:.6f}, shift {v.h100 - v.l40s:.6f}")


if __name__ == "__main__":
    main()
