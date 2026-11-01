"""Charts for the gpu-drift slides: what a call is, what a flip is, and how often calls flip between two GPUs.

Reads the recorded BRCA1 runs in a sibling evo2-cloud-operations checkout: the same 3,893 SNVs scored with Evo 2 7B
on an H100 (platform P2) and an L40S (platform P5), and the Youden threshold fitted on the H100 scores against the
Findlay assay labels (evidence/classification-flips/results.json). The score is the loss-of-function score,
-(variant - reference mean log-likelihood); higher means more damaging.

  uv run --with matplotlib --with pandas python gpu-drift/flip_charts.py [--brand gel]

Writes gpu_flip_step1, gpu_flip_step2 and gpu_flip_scatter (png, _gel for the GEL palette) into gpu-drift/out/.
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
OPS = HERE.parents[1] / "evo2-cloud-operations" / "evidence"
RUN = OPS / "runs" / "recorded-b75c67aa" / "platforms"
OUT = HERE / "out"
BRANDS = {
    "rewire": {"dark": "#16233A", "muted": "#5E6874", "line": "#7A8591", "surface": "#FDFCFA", "dot": "#2A6CB0",
               "flip": "#B4530F", "suffix": ""},
    "gel": {"dark": "#2B2F3B", "muted": "#60636C", "line": "#80828A", "surface": "#FFFFFF", "dot": "#005EB8",
            "flip": "#DF007D", "suffix": "_gel"},
}


def scores() -> tuple[pd.DataFrame, float]:
    h = pd.read_csv(RUN / "P2/attempt-1/w1_brca1_scores.tsv", sep="\t")
    l = pd.read_csv(RUN / "P5/attempt-2/w1_brca1_scores.tsv", sep="\t")
    d = h[["id", "delta"]].merge(l[["id", "delta"]], on="id", suffixes=("_h100", "_l40s"))
    d["h100"], d["l40s"] = -d["delta_h100"], -d["delta_l40s"]
    cut = json.loads((OPS / "classification-flips/results.json").read_text())["thresholds_full_set"]["youden"]
    d["flip"] = (d["h100"] >= cut) != (d["l40s"] >= cut)
    return d, cut


def number_line(c: dict, path: Path, v, cut: float, with_l40s: bool) -> None:
    fig, ax = plt.subplots(figsize=(14, 3.6))
    lo, hi = cut - 0.0006, cut + 0.0006
    ax.axhline(0, color=c["line"], lw=2, zorder=1)
    ax.axvline(cut, color=c["dark"], lw=2.5, ls="--", zorder=1)
    ax.text(cut + 0.00002, 0.82, "cut-off", ha="left", va="bottom", fontsize=20, color=c["dark"])
    ax.text(lo + 0.00003, -0.62, "functional", ha="left", fontsize=22, color=c["muted"])
    ax.text(hi - 0.00003, -0.62, "loss of function", ha="right", fontsize=22, color=c["muted"])
    ax.scatter([v.h100], [0], s=420, color=c["dot"], zorder=3)
    ax.text(v.h100, 0.32, "H100", ha="center", fontsize=20, color=c["dot"], fontweight="bold")
    if with_l40s:
        ax.scatter([v.l40s], [0], s=420, color=c["flip"], zorder=3)
        ax.text(v.l40s, 0.32, "L40S", ha="center", fontsize=20, color=c["flip"], fontweight="bold")
        ax.annotate("", xy=(v.l40s, -0.18), xytext=(v.h100, -0.18),
                    arrowprops=dict(arrowstyle="->", color=c["flip"], lw=2.5))
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.9, 1.25)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=c["surface"])
    plt.close(fig)


def scatter(c: dict, path: Path, d: pd.DataFrame, cut: float) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    lo, hi = cut - 0.0025, cut + 0.0035
    k = d[~d.flip]
    ax.scatter(k.h100, k.l40s, s=14, color=c["dot"], alpha=0.45, linewidths=0)
    f = d[d.flip]
    ax.scatter(f.h100, f.l40s, s=46, color=c["flip"], linewidths=0)
    ax.axvline(cut, color=c["dark"], lw=1.6, ls="--")
    ax.axhline(cut, color=c["dark"], lw=1.6, ls="--")
    ax.text(hi, cut + 0.00006, "cut-off", ha="right", va="bottom", fontsize=15, color=c["dark"])
    ax.text(lo + 0.0001, hi - 0.0001, f"{len(f)} of {len(d):,} calls flip", ha="left", va="top", fontsize=19,
            color=c["flip"], fontweight="bold")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Score on the H100", fontsize=17, color=c["dark"])
    ax.set_ylabel("Score on the L40S", fontsize=17, color=c["dark"])
    ax.tick_params(labelsize=12, colors=c["muted"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor(c["surface"])
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=c["surface"])
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--brand", default="rewire", choices=BRANDS)
    c = BRANDS[p.parse_args().brand]
    plt.rcParams.update({"font.family": "Arial"})
    OUT.mkdir(exist_ok=True)
    d, cut = scores()
    # The example: of the variants called loss of function on the H100 and functional on the L40S, the one whose
    # score moved most, so the shift is visible at this scale.
    flips = d[d.flip & (d.h100 >= cut)].assign(gap=lambda x: x.h100 - x.l40s).sort_values("gap", ascending=False)
    v = flips.iloc[0]
    number_line(c, OUT / f"gpu_flip_step1{c['suffix']}.png", v, cut, False)
    number_line(c, OUT / f"gpu_flip_step2{c['suffix']}.png", v, cut, True)
    scatter(c, OUT / f"gpu_flip_scatter{c['suffix']}.png", d, cut)
    print(f"cut-off {cut:.6f}; flips {int(d.flip.sum())} of {len(d)}; example {v.id}: H100 {v.h100:.6f}, "
          f"L40S {v.l40s:.6f}")


if __name__ == "__main__":
    main()
