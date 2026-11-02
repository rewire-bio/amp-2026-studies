"""Charts for the talk's slides, from the second run's results, in the talk's two brand palettes.

rd_charts.py draws the second run's own charts; this script draws the ones the slides use, sized for a 1920 x 1080
slide and coloured with the talk's palettes (rewire.bio, and Genomics England with files ending _gel). Dot and
interval throughout: a dot for the estimate, a line for its 95% interval (98.33% for the GPN family).

  v2_top1_scores      top 1 on the 500 for every single score and both baselines
  v2_by_group         top 1 on the 500 by diagnostic group: the best single lookup, Rule-based, Lookups combined
  v2_gains_avi        AVI added to each baseline, top 1: the 500 (locked) and the 200 (cross-validated)
  v2_gains_by_group   AVI added to each baseline by diagnostic group, top 1, on the 500
  v2_gains_gpn        the GPN family added to the lookups on the 500, top 1, with the 98.33% interval
  v2_first_vs_second  AVI added to the lookups: the first run's three results against this run's primary result
  v2_checks           AVI added to each baseline on the 500 under the main filters and each check

  cd rare-disease-v2 && uv run -q --with matplotlib python figures/slide_charts.py [--brand gel] [--out figures/out]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
RESULTS = HERE / "results"
FIRST = HERE.parent / "rare-disease" / "results"
BRANDS = {
    "rewire": {"grey": "#7A8591", "fm": "#2A6CB0", "evo": "#B4530F", "dark": "#16233A", "surface": "#FDFCFA",
               "font": "Arial", "suffix": ""},
    "gel": {"grey": "#8A8D93", "fm": "#005EB8", "evo": "#DF007D", "dark": "#2B2F3B", "surface": "#FFFFFF",
            "font": "Arial", "suffix": "_gel"},
}
GROUPS = [("canonical_splice", "Canonical splice", "SpliceAI"), ("nonsense", "Nonsense", "CADD"),
          ("missense", "Missense", "REVEL"), ("noncoding", "Non-coding", "SpliceAI")]
PRIMARY = "Lookups + AVI (locked) minus Lookups combined (locked)"
SECONDARY = "Rule + AVI minus Rule-based"


def load(name: str, folder: Path = RESULTS) -> dict:
    return json.loads((folder / name).read_text())


def pct(row: dict, key: str = "top_1", ci: str = "top_1_ci95") -> tuple[float, float, float]:
    return row[key] * 100, row[ci][0] * 100, row[ci][1] * 100


def num(x: float, signed: bool) -> str:
    """One decimal; never "-0.0" or "+0.0"."""
    x = round(x, 1) + 0.0
    return ("{:+.1f}" if signed and x != 0 else "{:.1f}").format(x)


def dots(ax, rows: list[tuple[str, tuple[float, float, float], str, bool]], c: dict, signed: bool) -> None:
    """One row per result, top to bottom: a dot for the estimate, a line for its interval, the value beside it."""
    for i, (label, (v, lo, hi), colour, filled) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=c[colour], lw=3, solid_capstyle="butt")
        ax.plot([v], [i], "o", ms=13, color=c[colour], mfc=c[colour] if filled else c["surface"], mew=2.5, zorder=3)
        ax.text(hi + 0.6, i, f"{num(v, signed)}  ({num(lo, signed)} to {num(hi, signed)})", va="center",
                fontsize=15, color=c["dark"])
    ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=17)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)


def gain_axes(ax, c: dict, lo: float, hi: float, label: str) -> None:
    ax.axvline(0, color=c["dark"], lw=1.4)
    for g in (-5, 5):
        ax.axvline(g, color=c["grey"], lw=1, ls=(0, (4, 4)))
    ax.set_xlim(lo, hi)
    ax.set_xlabel(label)


def save(fig, out: Path, name: str, c: dict) -> None:
    fig.tight_layout()
    fig.savefig(out / f"{name}{c['suffix']}.png", dpi=150)
    plt.close(fig)


def top1_scores(c: dict, out: Path) -> None:
    s = load("replication_main.json")["strata"]["all"]["scores"]
    rows = [("Lookups combined", s["Lookups combined (locked)"], "dark"), ("Rule-based", s["Rule-based"], "dark"),
            ("AVI alone", s["AVI"], "fm"), ("CADD", s["CADD"], "grey"), ("GPN-Star-V", s["GPN-Star-V"], "fm"),
            ("Evo 2 alone", s["Evo 2"], "evo"), ("SpliceAI", s["SpliceAI"], "grey"), ("GPN-MSA", s["GPN-MSA"], "fm"), ("REVEL", s["REVEL"], "grey"),
            ("AlphaMissense", s["AlphaMissense"], "grey")]
    fig, ax = plt.subplots(figsize=(14, 6.0))
    dots(ax, [(n, pct(r), col, True) for n, r, col in rows], c, signed=False)
    ax.set_xlim(0, 112)
    ax.set_xlabel("Spiked-in variant ranked first, 500 new cases (%)")
    save(fig, out, "v2_top1_scores", c)


def by_group(c: dict, out: Path) -> None:
    st = load("replication_main.json")["strata"]
    fig, axes = plt.subplots(1, len(GROUPS), figsize=(15, 4.2), sharey=True)
    for ax, (key, title, tool) in zip(axes, GROUPS):
        s = st[key]["scores"]
        rows = [("Its lookup alone", pct(s[tool]), "grey", True), ("Rule-based", pct(s["Rule-based"]), "dark", True),
                ("Lookups combined", pct(s["Lookups combined (locked)"]), "fm", True)]
        for i, (_, (v, lo, hi), colour, _) in enumerate(rows):
            ax.plot([lo, hi], [i, i], color=c[colour], lw=3, solid_capstyle="butt")
            ax.plot([v], [i], "o", ms=12, color=c[colour], zorder=3)
            label = f"{v:.0f} ({tool})" if i == 0 else f"{v:.0f}"
            ax.text(v, i - 0.32, label, ha="right" if v > 60 else "left", va="bottom", fontsize=15, color=c["dark"])
        ax.set_title(f"{title} ({st[key]['cases']})", fontsize=18, color=c["dark"], pad=10)
        ax.set_xlim(0, 105)
        ax.set_ylim(len(rows) - 0.4, -0.8)
        ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=16)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[0].set_xlabel("Ranked first (%)")
    save(fig, out, "v2_by_group", c)


def gains_avi(c: dict, out: Path) -> None:
    rep, snv = load("replication_main.json")["strata"]["all"]["contrasts"], load("snv_main.json")["strata"]["all"]["contrasts"]
    rows = [("Added to the lookups, 500 new cases", pct(rep[PRIMARY]), "fm", True),
            ("Added to the rule, 500 new cases", pct(rep[SECONDARY]), "dark", True),
            ("Added to the lookups, the 200 cases", pct(snv["Lookups + AVI minus Lookups combined"]), "fm", False),
            ("Added to the rule, the 200 cases", pct(snv[SECONDARY]), "dark", False)]
    fig, ax = plt.subplots(figsize=(14, 4.6))
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -8, 16, "Gain at top 1 from adding AVI (points); dashed lines at -5 and +5")
    save(fig, out, "v2_gains_avi", c)


def gains_by_group(c: dict, out: Path) -> None:
    st = load("replication_main.json")["strata"]
    rows = []
    for key, title in [("noncoding", "Non-coding"), ("noncoding: intronic", "Intronic"),
                       ("noncoding: splice_region", "Splice region"), ("missense", "Missense")]:
        n = st[key]["cases"]
        rows.append((f"{title} ({n}), to the lookups", pct(st[key]["contrasts"][PRIMARY]), "fm", True))
        rows.append((f"{title} ({n}), to the rule", pct(st[key]["contrasts"][SECONDARY]), "dark", True))
    fig, ax = plt.subplots(figsize=(14, 5.6))
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -10, 38, "Gain at top 1 from adding AVI, 500 new cases (points)")
    save(fig, out, "v2_gains_by_group", c)


def gains_gpn(c: dict, out: Path) -> None:
    con = load("replication_main.json")["strata"]["all"]["contrasts"]
    rows = [(f"{n} added to the lookups", pct(con[f"Lookups + {n} (locked) minus Lookups combined (locked)"],
                                             ci="top_1_ci98.33"), "fm", True)
            for n in ("GPN-MSA", "GPN-Star-V", "GPN-Star-M")]
    rows.append(("Evo 2 added to the lookups", pct(con["Lookups + Evo 2 (locked) minus Lookups combined (locked)"],
                                                   ci="top_1_ci98.75"), "evo", True))
    rows.append(("AVI added to the lookups (95%)", pct(con[PRIMARY]), "fm", True))
    fig, ax = plt.subplots(figsize=(14, 4.4))
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -8, 12, "Gain at top 1, 500 new cases (points)")
    save(fig, out, "v2_gains_added", c)


def indel_top1(c: dict, out: Path) -> None:
    st = load("indel_main.json")["strata"]["all"]
    s, con = st["scores"], st["contrasts"]
    rows = [("Lookups + Evo 2", s["Lookups + Evo 2"], "dark"), ("Lookups combined", s["Lookups combined"], "dark"),
            ("Evo 2 alone", s["Evo 2"], "evo"), ("Rule-based", s["Rule-based"], "dark"), ("SpliceAI", s["SpliceAI"], "grey")]
    fig, ax = plt.subplots(figsize=(14, 4.4))
    dots(ax, [(n, pct(r), col, True) for n, r, col in rows], c, signed=False)
    ax.set_xlim(0, 112)
    ax.set_xlabel("Spiked-in indel ranked first, 120 indel cases (%)")
    save(fig, out, "v2_indel_top1", c)


def first_vs_second(c: dict, out: Path) -> None:
    f200 = load("results.json", FIRST)["strata"]["all"]["added_to_lookups"]["Lookups + AlphaGenome AVI"]
    key = "Lookups + AlphaGenome AVI (locked)"
    f500 = load("results_replication.json", FIRST)["strata"]["all"]["added_to_locked_lookups"][key]
    fgrp = load("results_replication_gnomad0.001_grpmax.json", FIRST)["strata"]["all"]["added_to_locked_lookups"][key]
    v2 = load("replication_main.json")["strata"]["all"]["contrasts"][PRIMARY]
    rows = [("First run, first 200 cases", pct(f200), "grey", False),
            ("First run, 500 new cases", pct(f500), "grey", True),
            ("First run, 500 new, strict filter", pct(fgrp), "grey", True),
            ("Second run, 500 new cases", pct(v2), "fm", True)]
    fig, ax = plt.subplots(figsize=(14, 4.4))
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -8, 20, "Gain at top 1 from adding AVI to the lookups (points)")
    save(fig, out, "v2_first_vs_second", c)


def checks(c: dict, out: Path) -> None:
    rows = []
    for name, label in [("main", "Main filters"), ("bg1pct", "Looser filter (below 1%)"),
                        ("nomoi", "No inheritance filter"), ("all1kg", "All of 1000 Genomes left out")]:
        d = load(f"replication_{name}.json")
        con = d["strata"]["all"]["contrasts"]
        med = d["median_variants_per_case"]
        rows.append((f"{label} (median list {med:g}): to the lookups", pct(con[PRIMARY]), "fm", True))
        rows.append(("to the rule", pct(con[SECONDARY]), "dark", True))
    fig, ax = plt.subplots(figsize=(14, 5.8))
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -10, 16, "Gain at top 1 from adding AVI, 500 new cases (points)")
    save(fig, out, "v2_checks", c)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--brand", default="rewire", choices=BRANDS)
    p.add_argument("--out", type=Path, default=HERE / "figures" / "out")
    a = p.parse_args()
    c = BRANDS[a.brand]
    a.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": c["font"], "font.size": 22, "axes.edgecolor": "#C9CED4",
                         "axes.labelcolor": c["dark"], "xtick.color": c["dark"], "ytick.color": c["dark"],
                         "figure.facecolor": c["surface"], "axes.facecolor": c["surface"]})
    for chart in (top1_scores, by_group, gains_avi, gains_by_group, gains_gpn, indel_top1, checks):
        chart(c, a.out)
    print(f"wrote 7 charts to {a.out} ({a.brand})")


if __name__ == "__main__":
    main()
