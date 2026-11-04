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
# strong: the flipped call and the L40S; mid: the H100; rest: calls that do not flip. One chart style for the whole
# talk: Arial at slide sizes, the section's colours, a key row instead of a legend box, no grid, a light baseline.
BRANDS = {
    "rewire": {"dark": "#16233A", "muted": "#5E6874", "light": "#C9CED4", "surface": "#FDFCFA",
               "strong": "#2E8B7A", "mid": "#8CC9BC", "rest": "#C9CED4", "suffix": ""},
    "gel": {"dark": "#2B2F3B", "muted": "#60636C", "light": "#C3C3BF", "surface": "#FFFFFF",
            "strong": "#26913C", "mid": "#97E4A7", "rest": "#C3C3BF", "suffix": "_gel"},
}
# Charts are drawn at the size they are shown on a 1920 x 1080 slide: SLOT px per inch on screen, saved at twice
# that for sharpness. PX converts on-screen pixels to points, so text sizes below read as slide pixels.
SLOT, PX = 128, 72 / 128


def scores() -> tuple[pd.DataFrame, float]:
    h = pd.read_csv(RUN / "P2/attempt-1/w1_brca1_scores.tsv", sep="\t")
    l = pd.read_csv(RUN / "P5/attempt-2/w1_brca1_scores.tsv", sep="\t")
    d = h[["id", "delta"]].merge(l[["id", "delta"]], on="id", suffixes=("_h100", "_l40s"))
    d["h100"], d["l40s"] = -d["delta_h100"], -d["delta_l40s"]
    cut = json.loads((OPS / "classification-flips/results.json").read_text())["thresholds_full_set"]["youden"]
    d["flip"] = (d["h100"] >= cut) != (d["l40s"] >= cut)
    return d, cut


def number_line(c: dict, path: Path, v, cut: float, with_l40s: bool) -> None:
    fig, ax = plt.subplots(figsize=(1616 / SLOT, 430 / SLOT))
    lo, hi = cut - 0.0006, cut + 0.0006
    ax.axhline(0, color=c["light"], lw=4 * PX, zorder=1)
    ax.plot([cut, cut], [-0.55, 0.75], color=c["dark"], lw=3 * PX, ls=(0, (4, 3)), zorder=1)
    ax.text(cut, 0.82, "cut-off", ha="center", va="bottom", fontsize=28 * PX, color=c["dark"])
    ax.text(lo + 0.00002, -0.62, "functional", ha="left", va="top", fontsize=30 * PX, color=c["muted"])
    ax.text(hi - 0.00002, -0.62, "loss of function", ha="right", va="top", fontsize=30 * PX, color=c["muted"])
    ax.plot([v.h100], [0], "o", ms=34 * PX, color=c["mid"], mew=0, zorder=3)
    ax.text(v.h100, 0.3, "H100", ha="center", va="bottom", fontsize=30 * PX, color=c["dark"], fontweight="bold")
    if with_l40s:
        ax.plot([v.l40s], [0], "o", ms=34 * PX, color=c["strong"], mew=0, zorder=3)
        ax.text(v.l40s, 0.3, "L40S", ha="center", va="bottom", fontsize=30 * PX, color=c["dark"],
                fontweight="bold")
        ax.annotate("", xy=(v.l40s + 0.00003, -0.25), xytext=(v.h100 - 0.00003, -0.25),
                    arrowprops=dict(arrowstyle="-|>", color=c["strong"], lw=4 * PX, mutation_scale=30 * PX))
    ax.set_xlim(lo, hi)
    ax.set_ylim(-1.05, 1.25)
    ax.axis("off")
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.99)
    fig.savefig(path, dpi=SLOT * 2, facecolor=c["surface"])
    plt.close(fig)


def scatter(c: dict, path: Path, d: pd.DataFrame, cut: float) -> None:
    w, h = 740, 634
    fig, ax = plt.subplots(figsize=(w / SLOT, h / SLOT))
    lo, hi = cut - 0.0025, cut + 0.0035
    k, f = d[~d.flip], d[d.flip]
    ax.scatter(k.h100, k.l40s, s=(7 * PX) ** 2, color=c["rest"], linewidths=0, zorder=2)
    ax.scatter(f.h100, f.l40s, s=(12 * PX) ** 2, color=c["strong"], linewidths=0, zorder=3)
    ax.axvline(cut, color=c["dark"], lw=2.5 * PX, ls=(0, (4, 3)), zorder=1)
    ax.axhline(cut, color=c["dark"], lw=2.5 * PX, ls=(0, (4, 3)), zorder=1)
    ax.text(hi, cut + 0.00008, "cut-off", ha="right", va="bottom", fontsize=24 * PX, color=c["dark"])
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ticks = [0, 0.002, 0.004]
    ax.set_xticks(ticks, [f"{t:.3f}" for t in ticks])
    ax.set_yticks(ticks, [f"{t:.3f}" for t in ticks])
    ax.set_xlabel("Score on the H100", fontsize=26 * PX, color=c["dark"])
    ax.set_ylabel("Score on the L40S", fontsize=26 * PX, color=c["dark"])
    ax.tick_params(labelsize=22 * PX, colors=c["muted"], color=c["light"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(c["light"])
    ax.set_facecolor(c["surface"])
    key_h = 92
    fig.tight_layout(pad=0.3, rect=(0, 0, 1, 1 - key_h / h))
    # Key row: one round mark and name per group, two lines, from the left edge of the plotting area.
    x = ax.get_position().x0
    for i, (colour, name) in enumerate([("strong", f"Call flips: {len(f)} of {len(d):,}"),
                                        ("rest", "Same call on both GPUs")]):
        y = 1 - (22 + 40 * i) / h
        fig.text(x, y, "\u25cf", color=c[colour], fontsize=30 * PX, va="center", ha="left")
        fig.text(x + 34 / w, y, name, color=c["dark"], fontsize=26 * PX, va="center", ha="left",
                 fontweight="bold" if colour == "strong" else "normal")
    fig.savefig(path, dpi=SLOT * 2, facecolor=c["surface"])
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
