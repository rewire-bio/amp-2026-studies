"""The two EVEE results charts (evee-results and evee-results-2), drawn in one style so they agree.

  uv run --with matplotlib --with pandas --with pyarrow python evee/evee_results_cards.py \
      --evee ../evee-manuscript/artifacts

Same rows in the same order, same axis, three decimals, the value always right of its dot, one card per
consequence. Writes evee_results_step1.svg (missense) and evee_results_step2.svg (five other SNV classes) into
evee/out/; copy them into the matching slide folders of amp-2026-slides.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
INK, MUTED, CARD, EDGE = "#2B2F3B", "#60636C", "#E6EFF8", "#CCDFF1"
EVEE, ZERO, TOOL = "#07C5F5", "#9BE8FC", "#C3C3BF"
# (method in the table, label, colour); the same order on both slides
ROWS = [("Covpool probe", "EVEE probe", EVEE), ("AlphaMissense", "AlphaMissense", TOOL), ("REVEL", "REVEL", TOOL),
        ("CADD v1.7", "CADD v1.7", TOOL), ("GPN-MSA", "GPN-MSA", ZERO), ("Evo2 loss", "Evo 2 zero-shot", ZERO)]
LO, HI = 0.72, 1.04  # room for the value labels; ticks stop at 0.95


def draw(d, panels, rows, size, name, legend):
    fig = plt.figure(figsize=size, dpi=100)
    W, H = size[0] * 100, size[1] * 100
    left, gap, top, bottom = 230, 24, 70 if legend else 20, 60
    pw = (W - left - 20 - gap * (len(panels) - 1)) / len(panels)
    for k, (cons, title) in enumerate(panels):
        x0 = left + k * (pw + gap)
        ax = fig.add_axes([(x0 + 8) / W, bottom / H, (pw - 16) / W, (H - top - bottom - 44) / H])
        ax.set_zorder(2)
        fig.patches.append(FancyBboxPatch((x0 / W, (bottom - 40) / H), pw / W, (H - top - bottom + 40) / H,
                                          boxstyle="round,pad=0,rounding_size=0.01", transform=fig.transFigure,
                                          facecolor=CARD, edgecolor=EDGE, linewidth=2, zorder=0))
        fig.text((x0 + 14) / W, (H - top - 30) / H, title, fontsize=24, fontweight="bold", color=INK)
        sub = d[d.consequence == cons].set_index("method")
        for y, (method, label, colour) in enumerate(rows):
            if method in sub.index and sub.loc[method, "n_valid"] >= 100:
                v = sub.loc[method, "auroc"]
                ax.scatter([v], [y], s=260, color=colour, zorder=3)
                ax.annotate(f"{v:.3f}", (v, y), xytext=(12, 0), textcoords="offset points", va="center",
                            fontsize=20, fontweight="bold", color=INK)
        ax.set_xlim(LO, HI)
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_xticks([0.75, 0.85, 0.95])
        ax.tick_params(axis="x", colors=MUTED, labelsize=16, length=0)
        ax.set_yticks(range(len(rows)), [r[1] for r in rows] if k == 0 else [])
        ax.tick_params(axis="y", length=0, labelsize=22, colors=INK)
        if k == 0:
            for t in ax.get_yticklabels():
                t.set_fontweight("bold" if t.get_text() == "EVEE probe" else "normal")
        ax.set_facecolor("white")
        for s in ax.spines.values():
            s.set_visible(False)
        ax.grid(axis="x", color=EDGE, linewidth=1)
        ax.set_axisbelow(True)
    if legend:
        fig.legend(handles=[Line2D([], [], ls="", marker="o", markersize=14, color=c, label=l) for l, c in
                            [("EVEE probe", EVEE), ("Zero-shot DNA models", ZERO), ("Existing tools", TOOL)]],
                   loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False, fontsize=20)
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / f"{name}.svg", facecolor="white")
    fig.savefig(OUT / f"{name}.png", facecolor="white")
    plt.close(fig)
    print(f"{name}.svg")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--evee", type=Path, required=True)
    a = p.parse_args()
    plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none"})
    d = pd.read_feather(a.evee / "deconf_benchmark.feather")
    draw(d, [("Missense", "Missense")], ROWS, (14.4, 5.37), "evee_results_step1", True)
    others = [r for r in ROWS if r[0] not in ("AlphaMissense", "REVEL")]
    draw(d, [("Intronic", "Intronic"), ("Synonymous", "Synonymous"), ("Splice Acceptor", "Splice acceptor"),
             ("Splice Donor", "Splice donor"), ("Nonsense", "Nonsense")], others, (16.8, 4.8), "evee_results_step2",
         True)


if __name__ == "__main__":
    main()
