"""Charts for the case-benchmark slides, from rare-disease/results/results[_<study>].json.

  <prefix>_by_type.png   share of cases with the diagnostic variant in the top 5, by diagnostic variant type:
                         the best established lookup, the best foundation-model table, and Evo 2 7B
  <prefix>_combined.png  top 1 / 5 / 10 for the combined lookups, alone and with each foundation-model score added

  uv run --with matplotlib python rare-disease/figures/rd_charts.py                     # SNV study, prefix rd
  uv run --with matplotlib python rare-disease/figures/rd_charts.py --study indel       # prefix rd_indel
  ... --brand gel                                                            # GEL colours, files end _gel
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
BRANDS = {
    "rewire": {"grey": "#7A8591", "fm": "#2A6CB0", "evo": "#B4530F", "dark": "#16233A", "surface": "#FDFCFA",
               "font": "Arial", "suffix": ""},
    "gel": {"grey": "#8A8D93", "fm": "#005EB8", "evo": "#DF007D", "dark": "#2B2F3B", "surface": "#FFFFFF",
            "font": "Arial", "suffix": "_gel"},
}
LOOKUPS = ["CADD v1.7", "SpliceAI", "REVEL", "AlphaMissense"]
FM = ["GPN-MSA", "GPN-Star-V", "GPN-Star-M", "AlphaGenome AVI"]
LABELS = {"canonical_splice": "Canonical splice", "missense": "Missense", "nonsense": "Nonsense",
          "noncoding_or_synonymous": "Synonymous, UTR\nor intronic", "frameshift": "Frameshift",
          "inframe": "In-frame", "noncoding": "UTR or intronic"}
# The cancer study's "noncoding" group also takes synonymous SNVs (build_cases.py); the indel study's does not.
CANCER_LABELS = LABELS | {"noncoding": "Synonymous, UTR\nor intronic"}


def best(scores: dict, names: list[str], covered: dict) -> str | None:
    names = [n for n in names if n in scores and covered.get(n, 1) > 0]
    return max(names, key=lambda n: (scores[n]["top_5"], scores[n]["top_1"])) if names else None


def by_type(r: dict, c: dict, path: Path, run_label: str, labels: dict = LABELS) -> None:
    groups = [g for g in r["strata"] if g != "all"]
    fig, axes = plt.subplots(1, len(groups), figsize=(3.9 * len(groups), 5.4), sharey=True)
    for ax, key in zip(axes, groups):
        s = r["strata"][key]["scores"]
        bars = []
        for names, colour, label in ((LOOKUPS, c["grey"], "Best lookup"), (FM, c["fm"], "Best FM table")):
            b = best(s, names, r["coverage_diagnostic"])
            if b:
                ran = run_label and (run_label == "run" or key in ("frameshift", "inframe"))
                name = f"{b}, run" if b == "SpliceAI" and ran else b
                bars.append((f"{label}\n({name})", s[b]["top_5"], colour))
        if "Evo 2 7B" in s:
            bars.append(("Evo 2 7B\n(we ran it)", s["Evo 2 7B"]["top_5"], c["evo"]))
        for i, (_, v, colour) in enumerate(bars):
            ax.bar(i, v * 100, color=colour, width=0.66)
            ax.text(i, v * 100 + 2, f"{v * 100:.0f}%", ha="center", va="bottom", fontsize=14, color=c["dark"])
        ax.set_xticks(range(len(bars)), [b[0] for b in bars], fontsize=10.5)
        ax.set_title(f"{labels.get(key, key)} ({r['strata'][key]['cases']})", fontsize=15, color=c["dark"], pad=10)
        ax.set_ylim(0, 112)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Diagnostic variant in top 5 (%)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def combined(r: dict, c: dict, path: Path) -> None:
    s = r["strata"]["all"]["scores"]
    cov = r["coverage_diagnostic"]
    names = [("Lookups combined", c["grey"])] + [(n, c["fm"]) for n in s if n.startswith("Lookups + ") and "Evo 2" not in n
                                                and cov.get(n.replace("Lookups + ", ""), 1) > 0]
    if "Lookups + Evo 2" in s:
        names.append(("Lookups + Evo 2", c["evo"]))
    ks = [("top_1", "Top 1"), ("top_5", "Top 5"), ("top_10", "Top 10")]
    fig, ax = plt.subplots(figsize=(13, 5.2))
    width = 0.84 / len(names)
    alphas = [0.4, 0.55, 0.7, 1.0]
    fm_i = 0
    for j, (n, colour) in enumerate(names):
        alpha = 1.0
        if colour == c["fm"]:
            alpha, fm_i = alphas[fm_i % len(alphas)], fm_i + 1
        for i, (k, _) in enumerate(ks):
            v = s[n][k] * 100
            lo, hi = (x * 100 for x in s[n][f"{k}_ci95"])
            x = i + (j - (len(names) - 1) / 2) * width
            ax.bar(x, v, width=width * 0.9, color=colour, alpha=alpha, label=n.replace("Lookups + ", "+ ") if i == 0 else None)
            ax.errorbar(x, v, yerr=[[v - lo], [hi - v]], color="#45505F", capsize=2, lw=1)
            ax.text(x, hi + 1.5, f"{v:.0f}", ha="center", va="bottom", fontsize=10.5, color=c["dark"])
    ax.set_xticks(range(len(ks)), [k[1] for k in ks], fontsize=16)
    ax.set_ylabel(f"Cases, % (n = {r['strata']['all']['cases']})")
    ax.set_ylim(0, 115)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=len(names), loc="upper left", fontsize=11.5)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--study", default="snv")
    p.add_argument("--brand", default="rewire", choices=BRANDS)
    a = p.parse_args()
    c = BRANDS[a.brand]
    plt.rcParams.update({"font.family": c["font"], "font.size": 15, "axes.edgecolor": "#C9CED4",
                         "axes.labelcolor": c["dark"], "xtick.color": c["dark"], "ytick.color": c["dark"],
                         "figure.facecolor": c["surface"], "axes.facecolor": c["surface"]})
    r = json.loads((RESULTS / ("results.json" if a.study == "snv" else f"results_{a.study}.json")).read_text())
    prefix = "rd" if a.study == "snv" else f"rd_{a.study}"
    by_type(r, c, OUT / f"{prefix}_by_type{c['suffix']}.png", "run" if a.study == "indel" else "indels" if a.study == "cancer" else "",
            CANCER_LABELS if a.study == "cancer" else LABELS)
    combined(r, c, OUT / f"{prefix}_combined{c['suffix']}.png")
    print(f"wrote {prefix}_by_type{c['suffix']}.png, {prefix}_combined{c['suffix']}.png")


if __name__ == "__main__":
    main()
