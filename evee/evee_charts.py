"""EVEE charts, redrawn from EVEE's published figure data with a distinct colour per group and method.

  uv run --with matplotlib --with pandas --with pyarrow python evee/evee_charts.py \
      --evee ../evee-manuscript/artifacts [--brand gel]

Writes evee_confounding[_gel].svg and evee_results[_gel].svg into evee/out/, drawn for a 1234 x 460 px slot. Copy
them into the evee-confounding and evee-results slide folders of amp-2026-slides. --steps also writes
evee_results_step1 and _step2, with only the first one and three consequence rows filled in, for the stepped slides.
EVEE data: github.com/goodfire-ai/evee-manuscript, artifacts/*.feather (read only).
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
BRANDS = {
    "rewire": {"ink": "#16233A", "muted": "#5E6874", "grid": "#E3E6EA", "surface": "#FDFCFA", "suffix": "",
               "lof": "#B4530F", "mis": "#2A6CB0", "nc": "#2E8B7A",
               "evee": "#B4530F", "evo": "#E8A36B", "am": "#2A6CB0", "revel": "#6B4FA0", "cadd": "#2E8B7A"},
    "gel": {"ink": "#2B2F3B", "muted": "#5E6874", "grid": "#E4E4E1", "surface": "#F4F4F2", "suffix": "_gel",
            "lof": "#DF007D", "mis": "#005EB8", "nc": "#26913C",
            "evee": "#DF007D", "evo": "#F29CCB", "am": "#005EB8", "revel": "#26913C", "cadd": "#2B2F3B"},
}
W, H = 12.34, 4.6


def frame(ax, c, xlabel):
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(c["grid"])
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=c["grid"], linewidth=1)
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, labelpad=8, color=c["muted"])


def save(fig, c, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / f"{name}{c['suffix']}.svg", facecolor=c["surface"])
    plt.close(fig)
    print(f"{name}{c['suffix']}.svg")


def confounding(evee: Path, c: dict):
    d = pd.read_feather(evee / "dataset_characterization.feather")
    d = d[d.dataset.str.startswith("Naive")]
    groups = {"Splice Donor": "lof", "Nonsense": "lof", "Splice Acceptor": "lof", "Initiator Codon": "lof",
              "Missense": "mis", "Intronic": "nc", "5' UTR": "nc", "Synonymous": "nc"}
    d = d[d.consequence.isin(groups)].sort_values("pathogenic_rate")
    fig, ax = plt.subplots(figsize=(W, H))
    for y, (_, r) in enumerate(d.iterrows()):
        ax.barh(y, r.pathogenic_rate, height=0.62, color=c[groups[r.consequence]])
        ax.text(r.pathogenic_rate + 1, y, f"{r.pathogenic_rate:.1f}%   (n = {r['count']:,})", va="center",
                fontsize=13, color=c["ink"])
    ax.set_yticks(range(len(d)), d.consequence.str.replace("Codon", "codon"))
    ax.set_xlim(0, 128)
    ax.set_xticks([0, 25, 50, 75, 100])
    frame(ax, c, "Share labelled pathogenic in ClinVar SNVs (EVEE's full set, 834,000 variants)")
    ax.legend(handles=[Patch(color=c["lof"], label="Loss of function or canonical splice"),
                       Patch(color=c["mis"], label="Missense"),
                       Patch(color=c["nc"], label="Synonymous or non-coding")],
              loc="lower right", frameon=False, fontsize=13)
    fig.tight_layout()
    save(fig, c, "evee_confounding")


def results(evee: Path, c: dict, show: int = 6, name: str = "evee_results"):
    d = pd.read_feather(evee / "deconf_benchmark.feather")
    rows = ["Missense", "Synonymous", "Intronic", "Splice Donor", "Splice Acceptor", "Nonsense"]
    series = [("Covpool probe", "EVEE probe", c["evee"], "o", 150),
              ("Evo2 loss", "Evo 2 zero-shot", c["evo"], "o", 110),
              ("AlphaMissense", "AlphaMissense", c["am"], "D", 100),
              ("REVEL", "REVEL", c["revel"], "^", 120),
              ("CADD v1.7", "CADD v1.7", c["cadd"], "s", 100)]
    fig, ax = plt.subplots(figsize=(W, H))
    for r in range(len(rows)):
        if r % 2 == 0:
            ax.axhspan(r - 0.5, r + 0.5, color=c["grid"], alpha=0.35, zorder=0)
    for i, (method, label, colour, marker, size) in enumerate(series):
        sub = d[d.method == method].set_index("consequence")
        xs, ys = [], []
        for r, row in enumerate(rows):
            if r < show and row in sub.index and sub.loc[row, "n_valid"] >= 100:
                xs.append(sub.loc[row, "auroc"])
                ys.append(r + (i - 2) * 0.13)
        ax.scatter(xs, ys, s=size, marker=marker, color=colour, edgecolor=c["surface"], linewidth=1.2, zorder=3)
    ax.set_yticks(range(len(rows)), [r.capitalize() if r != "Splice Donor" else "Splice donor" for r in
                                     [x.replace("Splice Acceptor", "Splice acceptor") for x in rows]])
    ax.invert_yaxis()
    ax.set_xlim(0.7, 1.0)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    frame(ax, c, "AUROC within each consequence, deconfounded ClinVar set")
    handles = [Line2D([], [], linestyle="", marker=m, markersize=s ** 0.5, color=col, label=lab)
               for _, lab, col, m, s in series]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.17), frameon=False, fontsize=13,
              ncol=5, handletextpad=0.3, columnspacing=1.6)
    fig.tight_layout()
    save(fig, c, name)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--evee", type=Path, required=True)
    p.add_argument("--brand", default="rewire", choices=BRANDS)
    p.add_argument("--steps", action="store_true")
    a = p.parse_args()
    c = BRANDS[a.brand]
    plt.rcParams.update({"font.family": "Arial", "font.size": 14, "text.color": c["ink"],
                         "xtick.color": c["muted"], "ytick.color": c["ink"], "figure.facecolor": c["surface"],
                         "axes.facecolor": c["surface"], "svg.fonttype": "none"})
    confounding(a.evee, c)
    results(a.evee, c)
    for k, show in enumerate([1, 3] if a.steps else [], 1):
        results(a.evee, c, show, f"evee_results_step{k}")


if __name__ == "__main__":
    main()
