"""Charts for the talk's slides, from the second run's results, in the talk's two brand palettes.

rd_charts.py draws the second run's own charts; this script draws the ones the slides use, sized for a 1920 x 1080
slide and coloured with the talk's palettes (rewire.bio, and Genomics England with files ending _gel). Top 1 is a
ranked bar with its 95% interval as a thin line; a gain is a dot with its interval (98.33% for the GPN family, 98.75%
for Evo 2), over a band from -5 to +5 points, the pre-registered "no gain that matters".

  v2_top1_scores      top 1 on the 500 for every single score and both baselines
  v2_indel_top1       top 1 on the 120 indel cases, Evo 2 alone and added to the lookups
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
    "rewire": {"grey": "#7A8591", "light": "#C9CED4", "fm": "#2A6CB0", "evo": "#B4530F", "dark": "#16233A",
               "surface": "#FDFCFA", "strong": "#B4530F", "t1": "#F6E7DA", "t3": "#E8A36B",
               "font": "Arial", "suffix": ""},
    "gel": {"grey": "#60636C", "light": "#C3C3BF", "fm": "#005EB8", "evo": "#DF007D", "dark": "#2B2F3B",
            "surface": "#FFFFFF", "strong": "#DF007D", "t1": "#FFE2F2", "t3": "#FF8CCD",
            "font": "Arial", "suffix": "_gel"},
}
# Charts are drawn at the size they are shown on a 1920 x 1080 slide: SLOT px per inch on screen, saved at twice
# that for sharpness. PX converts on-screen pixels to points, so text sizes below read as slide pixels.
SLOT, PX = 128, 72 / 128
GROUPS = [("canonical_splice", "Canonical splice", "SpliceAI"), ("nonsense", "Nonsense", "CADD"),
          ("missense", "Missense", "REVEL"), ("noncoding", "Non-coding", "SpliceAI")]
PRIMARY = "Lookups + AVI (locked) minus Lookups combined (locked)"
SECONDARY = "Rule + AVI minus Rule-based"


def load(name: str, folder: Path = RESULTS) -> dict:
    return json.loads((folder / name).read_text())


def pct(row: dict, key: str = "top_1", ci: str = "top_1_ci95") -> tuple[float, float, float]:
    return row[key] * 100, row[ci][0] * 100, row[ci][1] * 100


def figure(w_px: int, h_px: int):
    return plt.subplots(figsize=(w_px / SLOT, h_px / SLOT))


def num(x: float, signed: bool) -> str:
    """One decimal, a true minus sign; never "-0.0" or "+0.0"."""
    x = round(x, 1) + 0.0
    return ("{:+.1f}" if signed and x != 0 else "{:.1f}").format(x).replace("-", "\u2212")


def dots(ax, rows: list[tuple[str, tuple[float, float, float], str, bool]], c: dict, signed: bool,
         pairs: list[str] | None = None) -> None:
    """One row per result, top to bottom: a dot for the estimate, a line for its interval, the value beside it.
    With pairs, rows come in twos (to the lookups, to the rule): one label per pair, centred on it, and the first
    pair's values say which is which."""
    for i, (label, (v, lo, hi), colour, filled) in enumerate(rows):
        tag = f"   {label}" if pairs and i < 2 else ""
        ax.plot([lo, hi], [i, i], color=c[colour], lw=4 * PX * 2, solid_capstyle="butt", zorder=3)
        ax.plot([v], [i], "o", ms=24 * PX, color=c[colour], mfc=c[colour] if filled else c["surface"],
                mew=4 * PX, zorder=4)
        ax.text(hi, i, f"   {num(v, signed)}  ({num(lo, signed)} to {num(hi, signed)}){tag}",
                va="center", fontsize=27 * PX, color=c["dark"])
    if pairs:
        ax.set_yticks([2 * k + 0.5 for k in range(len(pairs))], pairs, fontsize=28 * PX)
        for k in range(1, len(pairs)):
            ax.axhline(2 * k - 0.5, color=c["light"], lw=1.5 * PX * 2, zorder=0.5)
    else:
        ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=28 * PX)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(c["light"])
    ax.tick_params(axis="y", length=0, pad=12)
    ax.tick_params(axis="x", labelsize=26 * PX, color=c["light"])


def gain_axes(ax, c: dict, lo: float, hi: float, label: str, n: int) -> None:
    """The reading rule as a band from -5 to +5, labelled once above the first row, and a line at zero."""
    ax.axvspan(-5, 5, color=c["t1"], lw=0, zorder=0)
    ax.plot([0, 0], [n - 0.45, -0.55], color=c["dark"], lw=2 * PX * 2, zorder=1)
    ax.text(0, -0.9, "no gain that matters", ha="center", va="center", fontsize=26 * PX, color=c["dark"])
    ax.set_ylim(n - 0.45, -1.25)
    ax.set_xlim(lo, hi)
    ax.set_xticks([x for x in range(-10, 45, 5) if lo <= x <= hi])
    ax.set_xlabel(label, fontsize=26 * PX, color=c["grey"], labelpad=8)


def save(fig, out: Path, name: str, c: dict) -> None:
    fig.tight_layout(pad=0.3)
    fig.savefig(out / f"{name}{c['suffix']}.png", dpi=SLOT * 2, facecolor=c["surface"])
    plt.close(fig)


def bars(ax, rows: list[tuple[str, dict, str, bool]], c: dict, xlabel: str) -> None:
    """Ranked horizontal bars, best at the top: the value and its 95% interval (a thin line) at the bar end.
    rows: (label, result, fill colour key, foundation model). Foundation-model bars are outlined."""
    rows = sorted(rows, key=lambda r: -r[1]["top_1"])
    for i, (label, r, colour, fm) in enumerate(rows):
        v, lo, hi = pct(r)
        ax.barh(i, v, height=0.68, color=c[colour], edgecolor=c["dark"] if fm else c[colour],
                linewidth=3 * PX * 2 if fm else 0, zorder=2)
        ax.plot([lo, hi], [i, i], color=c["dark"], lw=2 * PX * 2, solid_capstyle="butt", zorder=3)
        ax.text(hi + 1.2, i, num(v, False), va="center", fontsize=28 * PX, color=c["dark"], fontweight="bold")
    ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=28 * PX)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.set_xlim(0, 108)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(c["light"])
    ax.tick_params(axis="y", length=0, pad=12)
    ax.tick_params(axis="x", labelsize=26 * PX, color=c["light"])
    ax.set_xlabel(xlabel, fontsize=26 * PX, color=c["grey"], labelpad=8)


def top1_scores(c: dict, out: Path) -> None:
    s = load("replication_main.json")["strata"]["all"]["scores"]
    rows = [("Lookups combined", s["Lookups combined (locked)"], "strong", False),
            ("Rule-based", s["Rule-based"], "strong", False),
            ("AVI alone", s["AVI"], "t3", True), ("CADD", s["CADD"], "t3", False),
            ("GPN-Star-V", s["GPN-Star-V"], "t3", True), ("Evo 2 alone", s["Evo 2"], "t3", True),
            ("SpliceAI", s["SpliceAI"], "t3", False), ("GPN-MSA", s["GPN-MSA"], "t3", True),
            ("REVEL", s["REVEL"], "t3", False), ("AlphaMissense", s["AlphaMissense"], "t3", False)]
    fig, ax = figure(1664, 600)
    bars(ax, rows, c, "Spiked-in variant ranked first, 500 new cases (%); thin line: 95% interval")
    ax.text(107, len(rows) - 1.5, "Outlined: foundation-model scores", ha="right", va="center",
            fontsize=26 * PX, color=c["dark"])
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
    rows = [("Added to the lookups, 500 new cases", pct(rep[PRIMARY]), "strong", True),
            ("Added to the rule, 500 new cases", pct(rep[SECONDARY]), "dark", True),
            ("Added to the lookups, the 200 cases", pct(snv["Lookups + AVI minus Lookups combined"]), "strong", False),
            ("Added to the rule, the 200 cases", pct(snv[SECONDARY]), "dark", False)]
    fig, ax = figure(1664, 547)
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -8, 16, "Gain at top 1 from adding AVI (points); open dots: the 200, cross-validated", len(rows))
    save(fig, out, "v2_gains_avi", c)


def gains_by_group(c: dict, out: Path) -> None:
    st = load("replication_main.json")["strata"]
    rows, pairs = [], []
    for key, title in [("noncoding", "Non-coding"), ("noncoding: intronic", "Intronic"),
                       ("noncoding: splice_region", "Splice region"), ("missense", "Missense")]:
        n = st[key]["cases"]
        pairs.append(f"{title}\n{n} cases")
        rows.append(("to the lookups", pct(st[key]["contrasts"][PRIMARY]), "strong", True))
        rows.append(("to the rule", pct(st[key]["contrasts"][SECONDARY]), "dark", True))
    fig, ax = figure(1664, 600)
    dots(ax, rows, c, signed=True, pairs=pairs)
    gain_axes(ax, c, -10, 42, "Gain at top 1 from adding AVI, 500 new cases (points)", len(rows))
    save(fig, out, "v2_gains_by_group", c)


def gains_gpn(c: dict, out: Path) -> None:
    con = load("replication_main.json")["strata"]["all"]["contrasts"]
    rows = [(f"{n} added to the lookups", pct(con[f"Lookups + {n} (locked) minus Lookups combined (locked)"],
                                             ci="top_1_ci98.33"), "strong", True)
            for n in ("GPN-MSA", "GPN-Star-V", "GPN-Star-M")]
    rows.append(("Evo 2 added to the lookups", pct(con["Lookups + Evo 2 (locked) minus Lookups combined (locked)"],
                                                   ci="top_1_ci98.75"), "strong", True))
    rows.append(("AVI added to the lookups (95%)", pct(con[PRIMARY]), "strong", True))
    fig, ax = figure(1664, 523)
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -8, 12, "Gain at top 1, 500 new cases (points); intervals 98.33% for GPN, 98.75% for Evo 2",
              len(rows))
    save(fig, out, "v2_gains_added", c)


def indel_top1(c: dict, out: Path) -> None:
    s = load("indel_main.json")["strata"]["all"]["scores"]
    rows = [("Lookups + Evo 2", s["Lookups + Evo 2"], "strong", False),
            ("Lookups combined", s["Lookups combined"], "strong", False),
            ("Evo 2 alone", s["Evo 2"], "t3", True), ("Rule-based", s["Rule-based"], "strong", False),
            ("SpliceAI", s["SpliceAI"], "t3", False)]
    fig, ax = figure(1728, 560)
    bars(ax, rows, c, "Spiked-in indel ranked first, 120 indel cases (%); thin line: 95% interval")
    ax.text(107, len(rows) - 1, "Outlined: foundation-model score", ha="right", va="center",
            fontsize=26 * PX, color=c["dark"])
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
    fig, ax = figure(1664, 523)
    dots(ax, rows, c, signed=True)
    gain_axes(ax, c, -8, 20, "Gain at top 1 from adding AVI to the lookups (points)", len(rows))
    save(fig, out, "v2_first_vs_second", c)


def checks(c: dict, out: Path) -> None:
    rows, pairs = [], []
    for name, label in [("main", "Main filters"), ("bg1pct", "Looser filter (below 1%)"),
                        ("nomoi", "No inheritance filter"), ("all1kg", "All of 1000 Genomes left out")]:
        d = load(f"replication_{name}.json")
        con = d["strata"]["all"]["contrasts"]
        med = d["median_variants_per_case"]
        pairs.append(f"{label}\nmedian list {med:g}")
        rows.append(("to the lookups", pct(con[PRIMARY]), "strong", True))
        rows.append(("to the rule", pct(con[SECONDARY]), "dark", True))
    fig, ax = figure(1664, 500)
    dots(ax, rows, c, signed=True, pairs=pairs)
    gain_axes(ax, c, -10, 24, "Gain at top 1 from adding AVI, 500 new cases (points)", len(rows))
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
