"""Charts for the case-benchmark slides, from rare-disease/results/results[_<study>].json.

  <prefix>_by_type.png       share of cases with the diagnostic variant ranked first, all cases and by diagnostic
                             variant type, for the same four series in every panel: the lookups combined, AVI alone,
                             the lookups + AVI and the lookups + Evo 2 (SNV and cancer studies)
  <prefix>_combined.png      top 1 / 5 / 10 for the combined lookups, alone and with AVI, GPN-MSA or Evo 2 added
  rd_gains.png               gain at top 1 from adding each score, first 200 cases against 500 new (model locked)
  rd_stricter.png            gain at top 1 from adding AVI as the tests get stricter (replication, frequency filters)

--steps also writes <name>_step<k>.png for the slides that build a chart up one part at a time: the later panels
(by type) or groups (combined) are left blank, so every step has the same size and layout.

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
LABELS = {"canonical_splice": "Canonical splice", "missense": "Missense", "nonsense": "Nonsense",
          "noncoding_or_synonymous": "Synonymous, UTR\nor intronic", "frameshift": "Frameshift",
          "inframe": "In-frame", "noncoding": "UTR or intronic"}
# The cancer study's "noncoding" group also takes synonymous SNVs (build_cases.py); the indel study's does not.
CANCER_LABELS = LABELS | {"noncoding": "Synonymous, UTR\nor intronic"}


SERIES = [("Lookups combined", "grey", 1.0, "Lookups combined"), ("AlphaGenome AVI", "fm", 0.45, "AVI alone"),
          ("Lookups + AlphaGenome AVI", "fm", 1.0, "Lookups + AVI"), ("Lookups + Evo 2", "evo", 1.0, "Lookups + Evo 2")]
ORDER = {"snv": ["all", "canonical_splice", "nonsense", "missense", "noncoding_or_synonymous"],
         "cancer": ["all", "missense", "canonical_splice", "nonsense", "frameshift", "noncoding"]}
STEPS = {"snv": [1, 3, 4, 5], "cancer": [4, 5, 6]}


def fixed_by_type(r: dict, c: dict, path: Path, groups: list[str], labels: dict, show: int) -> None:
    """Top 1 by type, the same series in every panel; panels after the first `show` are left blank."""
    series = [x for x in SERIES if x[0] in r["strata"]["all"]["scores"]]
    fig, axes = plt.subplots(1, len(groups), figsize=(3.6 * len(groups), 5.6), sharey=True)
    for k, (ax, key) in enumerate(zip(axes, groups)):
        if k >= show:
            ax.set_visible(False)
            continue
        s = r["strata"][key]["scores"]
        for i, (name, colour, alpha, _) in enumerate(series):
            v = s[name]["top_1"] * 100
            ax.bar(i, v, color=c[colour], alpha=alpha, width=0.72)
            ax.text(i, v + 2, f"{v:.0f}", ha="center", va="bottom", fontsize=17, color=c["dark"])
        ax.set_xticks([])
        title = "All cases" if key == "all" else labels.get(key, key)
        ax.set_title(f"{title} ({r['strata'][key]['cases']})", fontsize=19, color=c["dark"], pad=12,
                     fontweight="bold" if key == "all" else "normal")
        ax.set_ylim(0, 112)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Diagnostic variant ranked first (%)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c[col], alpha=a) for _, col, a, _ in series]
    fig.legend(handles, [x[3] for x in series], loc="lower center", ncol=len(series), frameon=False, fontsize=18,
               bbox_to_anchor=(0.5, 0))
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def combined_steps(r: dict, c: dict, path: Path, show: int) -> None:
    """Top 1 / 5 / 10: the lookups alone and with AVI, GPN-MSA or Evo 2; groups after the first `show` blank."""
    s = r["strata"]["all"]["scores"]
    names = [(n, col, a) for n, col, a in (("Lookups combined", "grey", 1.0), ("Lookups + AlphaGenome AVI", "fm", 1.0),
                                          ("Lookups + GPN-MSA", "fm", 0.45), ("Lookups + Evo 2", "evo", 1.0))
             if n in s and not (n == "Lookups + GPN-MSA" and r["coverage_diagnostic"].get("GPN-MSA", 1) == 0)]
    ks = [("top_1", "Top 1"), ("top_5", "Top 5"), ("top_10", "Top 10")]
    fig, ax = plt.subplots(figsize=(14, 5.2))
    width = 0.8 / len(names)
    for j, (n, colour, alpha) in enumerate(names):
        for i, (k, _) in enumerate(ks[:show]):
            v = s[n][k] * 100
            lo, hi = (x * 100 for x in s[n][f"{k}_ci95"])
            x = i + (j - (len(names) - 1) / 2) * width
            ax.bar(x, v, width=width * 0.9, color=c[colour], alpha=alpha, label=n.replace("Lookups + AlphaGenome AVI", "+ AVI").replace("Lookups + ", "+ ") if i == 0 else None)
            ax.errorbar(x, v, yerr=[[v - lo], [hi - v]], color="#45505F", capsize=2, lw=1)
            ax.text(x, hi + 1.5, f"{v:.0f}", ha="center", va="bottom", fontsize=16, color=c["dark"])
    ax.set_xticks(range(len(ks)), [k[1] if i < show else "" for i, k in enumerate(ks)], fontsize=21)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel(f"Cases, % (n = {r['strata']['all']['cases']})")
    ax.set_ylim(0, 115)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=len(names), loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=16, handlelength=1.2,
              columnspacing=1.4)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def gain_bars(ax, rows: list[tuple[str, float, list[float], str, float]], c: dict) -> None:
    for i, (label, v, ci, colour, alpha) in enumerate(rows):
        v, lo, hi = v * 100, ci[0] * 100, ci[1] * 100
        ax.bar(i, v, color=c[colour], alpha=alpha, width=0.66)
        ax.errorbar(i, v, yerr=[[v - lo], [hi - v]], color="#45505F", capsize=4, lw=1.4)
        ax.text(i, hi + 0.8, f"{v:+.1f}", ha="center", va="bottom", fontsize=18, color=c["dark"])
    ax.axhline(0, color="#45505F", lw=1.2)
    ax.set_xticks(range(len(rows)), [x[0] for x in rows], fontsize=16)
    ax.set_ylabel("Gain at top 1 (points)")
    ax.spines[["top", "right"]].set_visible(False)


def gains(c: dict, path: Path) -> None:
    """Added to the lookups, gain at top 1: first 200 cases against 500 new cases with the model locked."""
    first = json.loads((RESULTS / "results.json").read_text())["strata"]["all"]["added_to_lookups"]
    rep = json.loads((RESULTS / "results_replication.json").read_text())["strata"]["all"]["added_to_locked_lookups"]
    rows = []
    for name, short, colour in (("AlphaGenome AVI", "AVI", "fm"), ("GPN-MSA", "GPN-MSA", "fm"), ("Evo 2", "Evo 2", "evo")):
        a = first[f"Lookups + {name}"]
        rows.append((f"{short}\n200 cases", a["top_1"], a["top_1_ci95"], colour, 0.45))
        b = rep.get(f"Lookups + {name} (locked)")
        if b:
            rows.append((f"{short}\n500 new", b["top_1"], b["top_1_ci95"], colour, 1.0))
    fig, ax = plt.subplots(figsize=(14, 5.2))
    gain_bars(ax, rows, c)
    ax.set_ylim(-6, 20)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def stricter(c: dict, path: Path) -> None:
    """AVI added to the lookups, gain at top 1, as the tests get stricter."""
    runs = [("200 cases", "results.json", "added_to_lookups", ""),
            ("500 new", "results_replication.json", "added_to_locked_lookups", " (locked)"),
            ("500 new,\ngnomAD AF\nbelow 0.1%", "results_replication_gnomad0.001.json", "added_to_locked_lookups", " (locked)"),
            ("500 new,\nstrict filter:\ngrpmax below 0.1%", "results_replication_gnomad0.001_grpmax.json",
             "added_to_locked_lookups", " (locked)")]
    rows = []
    for i, (label, f, key, suffix) in enumerate(runs):
        a = json.loads((RESULTS / f).read_text())["strata"]["all"][key][f"Lookups + AlphaGenome AVI{suffix}"]
        rows.append((label, a["top_1"], a["top_1_ci95"], "fm", 0.45 if i == len(runs) - 1 else 1.0))
    fig, ax = plt.subplots(figsize=(14, 5.6))
    gain_bars(ax, rows, c)
    ax.set_ylim(-6, 20)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--study", default="snv")
    p.add_argument("--brand", default="rewire", choices=BRANDS)
    p.add_argument("--steps", action="store_true", help="also write the step images for stepped slides")
    a = p.parse_args()
    c = BRANDS[a.brand]
    plt.rcParams.update({"font.family": c["font"], "font.size": 18, "axes.edgecolor": "#C9CED4",
                         "axes.labelcolor": c["dark"], "xtick.color": c["dark"], "ytick.color": c["dark"],
                         "figure.facecolor": c["surface"], "axes.facecolor": c["surface"]})
    p_steps = a.steps
    r = json.loads((RESULTS / ("results.json" if a.study == "snv" else f"results_{a.study}.json")).read_text())
    prefix = "rd" if a.study == "snv" else f"rd_{a.study}"
    written = []
    if a.study in ORDER:
        labels = CANCER_LABELS if a.study == "cancer" else LABELS
        groups = ORDER[a.study]
        fixed_by_type(r, c, OUT / f"{prefix}_by_type{c['suffix']}.png", groups, labels, len(groups))
        written.append(f"{prefix}_by_type")
        for k, show in enumerate(STEPS[a.study] if p_steps else [], 1):
            fixed_by_type(r, c, OUT / f"{prefix}_by_type_step{k}{c['suffix']}.png", groups, labels, show)
            written.append(f"{prefix}_by_type_step{k}")
    if a.study != "cancer":
        combined_steps(r, c, OUT / f"{prefix}_combined{c['suffix']}.png", 3)
        written.append(f"{prefix}_combined")
        for k, show in enumerate([1] if p_steps else [], 1):
            combined_steps(r, c, OUT / f"{prefix}_combined_step{k}{c['suffix']}.png", show)
            written.append(f"{prefix}_combined_step{k}")
    if a.study == "snv":
        gains(c, OUT / f"rd_gains{c['suffix']}.png")
        stricter(c, OUT / f"rd_stricter{c['suffix']}.png")
        written += ["rd_gains", "rd_stricter"]
    print("wrote " + ", ".join(f"{w}{c['suffix']}.png" for w in written))


if __name__ == "__main__":
    main()
