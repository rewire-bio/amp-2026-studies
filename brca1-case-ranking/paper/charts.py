"""Charts for the BRCA1 case-ranking paper and companion, drawn from ../results.json only.

  uv run --with matplotlib python brca1-case-ranking/paper/charts.py

Writes figures/<name>.pdf (for the paper) next to this file and ../figures/<name>.png (for the companion page).
Every value comes from results.json; nothing is recomputed.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
RESULTS = json.loads((STUDY / "results.json").read_text())
PDF_DIR, PNG_DIR = HERE / "figures", STUDY / "figures"

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#8a8985"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"

plt.rcParams.update({"font.family": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 9,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "figure.dpi": 150})

ORDER = ["REVEL", "SpliceAI", "Evo 2 7B", "CADD v1.7",
         "Average percentile CADD + SpliceAI", "Average percentile CADD + SpliceAI + Evo 2",
         "CADD + SpliceAI", "CADD + Evo 2", "CADD + SpliceAI + Evo 2",
         "Type-aware CADD + SpliceAI", "Type-aware CADD + SpliceAI + Evo 2"]
LABEL = {"Average percentile CADD + SpliceAI": "Avg percentile C + S",
         "Average percentile CADD + SpliceAI + Evo 2": "Avg percentile C + S + Evo 2",
         "CADD + SpliceAI": "Model C + S", "CADD + Evo 2": "Model C + Evo 2",
         "CADD + SpliceAI + Evo 2": "Model C + S + Evo 2",
         "Type-aware CADD + SpliceAI": "Type-aware C + S",
         "Type-aware CADD + SpliceAI + Evo 2": "Type-aware C + S + Evo 2"}


def pct(x: float) -> float:
    return 100 * x


def save(fig, name: str) -> None:
    PDF_DIR.mkdir(exist_ok=True)
    PNG_DIR.mkdir(exist_ok=True)
    fig.savefig(PDF_DIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(PNG_DIR / f"{name}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("wrote", name)


def colour(name: str) -> str:
    return ORANGE if "Evo 2" in name else BLUE


def dot_panel(ax, rows, metric: str, title: str) -> None:
    """rows: list of (label, colour, result dict). Draw point and 95% interval per row, top row first."""
    n = len(rows)
    for i, (label, col, r) in enumerate(rows):
        y = n - 1 - i
        lo, hi = (pct(v) for v in r[f"{metric}_ci95"])
        ax.plot([lo, hi], [y, y], color=col, lw=2, solid_capstyle="round")
        ax.plot(pct(r[metric]), y, "o", color=col, ms=6, mec="white", mew=1)
        ax.text(hi + 1.5, y, f"{pct(r[metric]):.1f}", va="center", fontsize=8, color=MUTED)
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for lab, _, _ in rows][::-1])
    ax.set_xlim(0, 110)
    ax.set_xticks(range(0, 101, 20))
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=9.5)
    ax.set_xlabel("Cases (%), 95% interval")


def scores_k100() -> None:
    s = RESULTS["cases"]["k100"]["all"]["scores"]
    rows = [(LABEL.get(n, n), colour(n), s[n]) for n in ORDER]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    dot_panel(axes[0], rows, "ranked_first", "Ranked first")
    dot_panel(axes[1], rows, "in_top_5", "In the top 5")
    for ax in axes:
        ax.axhline(len(ORDER) - 4.5, color=GRID, lw=1)
    fig.suptitle("Where the loss-of-function variant lands among 101 candidates (8,230 cases)\n"
                 "Blue: no Evo 2. Orange: includes Evo 2. C = CADD v1.7, S = SpliceAI.",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    save(fig, "scores_k100")


PAIRS = [("CADD + Evo 2 vs CADD v1.7", "CADD alone"),
         ("Average percentile CADD + SpliceAI + Evo 2 vs Average percentile CADD + SpliceAI", "Avg percentile C + S"),
         ("CADD + SpliceAI + Evo 2 vs CADD + SpliceAI", "Model C + S"),
         ("Type-aware CADD + SpliceAI + Evo 2 vs Type-aware CADD + SpliceAI", "Type-aware C + S")]


def gains_k100() -> None:
    k = RESULTS["cases"]["k100"]
    strata = [("all", "All cases (8,230)"), ("coding", "Coding (5,840)"), ("noncoding", "Noncoding (2,390)")]
    fig, axes = plt.subplots(1, 2, figsize=(9, 5.2), sharey=True)
    for ax, metric, title in zip(axes, ["ranked_first", "in_top_5"], ["Ranked first", "In the top 5"]):
        labels, y = [], 0
        ticks = []
        for si, (stratum, slabel) in enumerate(strata):
            for key, blabel in PAIRS:
                d = k[stratum]["paired"][key][metric]
                lo, hi = (100 * v for v in d["ci95"])
                est = 100 * d["difference"]
                col = ORANGE if lo > 0 else BLUE if hi < 0 else GREY
                ax.plot([lo, hi], [-y, -y], color=col, lw=2, solid_capstyle="round")
                ax.plot(est, -y, "o", color=col, ms=6, mec="white", mew=1)
                ax.text(17, -y, f"{est:+.1f} ({lo:+.1f} to {hi:+.1f})", va="center", fontsize=7.5, color=MUTED)
                ticks.append(-y)
                labels.append(f"{slabel}: {blabel}" if blabel == PAIRS[0][1] else blabel)
                y += 1
            y += 0.6
        ax.axvline(0, color=MUTED, lw=1)
        ax.set_xlim(-12, 36)
        ax.set_xticks(range(-10, 16, 5))
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_title(title, loc="left", fontsize=9.5)
        ax.set_xlabel("Change from adding Evo 2, points (95% interval)")
    axes[0].set_yticks(ticks)
    axes[0].set_yticklabels(labels)
    fig.suptitle("Gain from adding Evo 2 to each baseline, 100 candidates per case\n"
                 "Orange: interval above 0. Blue: interval below 0. Grey: interval includes 0.",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    save(fig, "gains_k100")


def types_k100() -> None:
    k = RESULTS["cases"]["k100"]
    names = ["SpliceAI", "Evo 2 7B", "CADD v1.7", "CADD + Evo 2",
             "Type-aware CADD + SpliceAI", "Type-aware CADD + SpliceAI + Evo 2"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.0), sharey=True)
    for ax, stratum, title in zip(axes, ["coding", "noncoding"],
                                  [f"Coding diagnostic variant ({k['coding']['cases']:,} cases)",
                                   f"Noncoding diagnostic variant ({k['noncoding']['cases']:,} cases)"]):
        rows = [(LABEL.get(n, n), colour(n), k[stratum]["scores"][n]) for n in names]
        dot_panel(ax, rows, "in_top_5", title)
    fig.suptitle("In the top 5 by variant type, 100 candidates per case. Blue: no Evo 2. Orange: includes Evo 2.",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    save(fig, "types_k100")


if __name__ == "__main__":
    scores_k100()
    gains_k100()
    types_k100()
