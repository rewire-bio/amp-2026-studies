"""Charts for the talk's slides, from the second run's results, in the talk's two brand palettes.

rd_charts.py draws the second run's own charts; this script draws the ones the slides use, sized for a 1920 x 1080
slide and coloured with the talk's palettes (rewire.bio, and Genomics England with files ending _gel). One chart style
for the whole talk: horizontal bars, best at the top; colour shows the group (the section's strong colour for the
highlighted group, its mid tint for the second, warm grey for the rest), named in a key row above the chart; the value
in bold at the bar end; an interval as a thin line. A gain is a round dot with its interval (98.33% for the GPN family,
98.75% for Evo 2), over a band from -5 to +5 percentage points: the smallest gain worth having, set in advance, so a
gain inside it is too small to count.

  v2_top1_scores      top 1 on the 500 for every single score and both baselines
  v2_indel_top1       top 1 on the 120 indel cases, Evo 2 alone and added to the combined-score baseline
  v2_by_group         top 1 on the 500 by diagnostic group: the group's existing score alone and both baselines
  v2_gains_avi        AVI added to each baseline, top 1: the 500 (locked) and the 200 (cross-validated)
  v2_gains_by_group   AVI added to each baseline by diagnostic group, top 1, on the 500
  v2_gains_added      the GPN family and Evo 2 added to the combined-score baseline on the 500, top 1
  v2_first_vs_second  AVI added to the combined-score baseline: the first run's three results against this run's
  v2_checks           AVI added to each baseline on the 500 under the main filters and each check
  af_sweep            top 1 on the 500 against the grpmax threshold (frequency sweep, exploratory), every single
                      score and both baselines; PNG and SVG
  af_sweep_gain       AVI added to each baseline at each threshold of the sweep; PNG and SVG

The two sweep charts are drawn only when results/replication_sweep.json exists (case_benchmark.py --analysis sweep).
RDV2_RESULTS points the script at another results folder, for tests.

  cd rare-disease-v2 && uv run -q --with matplotlib python figures/slide_charts.py [--brand gel] [--out figures/out]
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
RESULTS = Path(os.environ.get("RDV2_RESULTS", HERE / "results"))
FIRST = HERE.parent / "rare-disease" / "results"
# strong: the highlighted group; mid: the second group; rest: everything else; band: the too-small-to-count band
BRANDS = {
    "rewire": {"grey": "#5E6874", "light": "#C9CED4", "dark": "#16233A", "surface": "#FDFCFA",
               "strong": "#B4530F", "mid": "#E8A36B", "rest": "#C9CED4", "band": "#F6E7DA", "suffix": ""},
    "gel": {"grey": "#60636C", "light": "#DBDCD9", "dark": "#2B2F3B", "surface": "#FFFFFF",
            "strong": "#DF007D", "mid": "#FF8CCD", "rest": "#C3C3BF", "band": "#FFE2F2", "suffix": "_gel"},
}
# Charts are drawn at the size they are shown on a 1920 x 1080 slide: SLOT px per inch on screen, saved at twice
# that for sharpness. PX converts on-screen pixels to points, so text sizes below read as slide pixels.
SLOT, PX = 128, 72 / 128
LABEL, VALUE, TICK, KEY = 30 * PX, 30 * PX, 26 * PX, 28 * PX
KEY_H = 60  # slide pixels kept above the chart for the key row
GROUPS = [("canonical_splice", "Canonical splice", "SpliceAI"), ("nonsense", "Nonsense", "CADD"),
          ("missense", "Missense", "REVEL"), ("noncoding", "Non-coding", "SpliceAI")]
PRIMARY = "Lookups + AVI (locked) minus Lookups combined (locked)"
SECONDARY = "Rule + AVI minus Rule-based"
TO_COMBINED, TO_RULE = "Added to the combined-score baseline", "Added to the rule baseline"


def load(name: str, folder: Path = RESULTS) -> dict:
    return json.loads((folder / name).read_text())


def pct(row: dict, key: str = "top_1", ci: str = "top_1_ci95") -> tuple[float, float, float]:
    return row[key] * 100, row[ci][0] * 100, row[ci][1] * 100


def figure(w_px: int, h_px: int, ncols: int = 1):
    return plt.subplots(1, ncols, figsize=(w_px / SLOT, h_px / SLOT), sharey=ncols > 1)


def num(x: float, signed: bool) -> str:
    """One decimal, a true minus sign; never "-0.0" or "+0.0"."""
    x = round(x, 1) + 0.0
    return ("{:+.1f}" if signed and x != 0 else "{:.1f}").format(x).replace("-", "−")


def style(ax, c: dict) -> None:
    """White ground, one light baseline, no grid, no border."""
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(c["light"])
    ax.spines["bottom"].set_linewidth(2 * PX)
    ax.tick_params(axis="y", length=0, pad=14, labelsize=LABEL, labelcolor=c["dark"])
    ax.tick_params(axis="x", labelsize=TICK, colors=c["grey"], color=c["light"])


def save(fig, out: Path, name: str, c: dict, key: list[tuple[str, str]], marker: str = "■", svg: bool = False) -> None:
    """Lay the chart out below a key row, then draw the key: one coloured mark and name per group, from the left
    edge of the plotting area. svg: also save an SVG."""
    h_px = fig.get_figheight() * SLOT
    fig.tight_layout(pad=0.3, rect=(0, 0, 1, 1 - KEY_H / h_px))
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    x, y = min(ax.get_position().x0 for ax in fig.axes), 1 - KEY_H / 2 / h_px
    for colour, group in key:
        t = fig.text(x, y, marker, color=c[colour], fontsize=KEY * 1.15, va="center", ha="left")
        x = t.get_window_extent(r).transformed(fig.transFigure.inverted()).x1 + 8 / (fig.get_figwidth() * SLOT)
        t = fig.text(x, y, group, color=c["dark"], fontsize=KEY, va="center", ha="left")
        x = t.get_window_extent(r).transformed(fig.transFigure.inverted()).x1 + 44 / (fig.get_figwidth() * SLOT)
    fig.savefig(out / f"{name}{c['suffix']}.png", dpi=SLOT * 2, facecolor=c["surface"])
    if svg:
        fig.savefig(out / f"{name}{c['suffix']}.svg", facecolor=c["surface"])
    plt.close(fig)


def bars(ax, rows: list[tuple[str, dict, str]], c: dict, xlabel: str) -> None:
    """Ranked horizontal bars, best at the top: the value in bold and its 95% interval (a thin grey line) at the bar
    end. rows: (label, result, colour key)."""
    rows = sorted(rows, key=lambda r: -r[1]["top_1"])
    for i, (label, r, colour) in enumerate(rows):
        v, lo, hi = pct(r)
        ax.barh(i, v, height=0.66, color=c[colour], linewidth=0, zorder=2)
        ax.plot([lo, hi], [i, i], color=c["grey"], lw=3 * PX, solid_capstyle="round", zorder=3)
        ax.text(hi + 1.2, i, num(v, False), va="center", fontsize=VALUE, color=c["dark"], fontweight="bold")
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.set_ylim(len(rows) - 0.5, -0.55)
    ax.set_xlim(0, 108)
    ax.set_xticks([0, 25, 50, 75, 100])
    style(ax, c)
    ax.set_xlabel(xlabel, fontsize=TICK, color=c["grey"], labelpad=8)


def dots(ax, rows: list[tuple[str, tuple[float, float, float], str]], c: dict, pairs: list[str] | None = None) -> list[float]:
    """One row per result, top to bottom: a round dot for the gain, a thin line in its colour for the interval, the
    value beside it. With pairs, rows come in twos (to the combined-score baseline, to the rule baseline) under one
    label per pair, with a gap between pairs."""
    step = 2.7 if pairs else 1
    ys = [(i // 2) * step + i % 2 if pairs else i for i in range(len(rows))]
    for y, (_, (v, lo, hi), colour) in zip(ys, rows):
        ax.plot([lo, hi], [y, y], color=c[colour], lw=4 * PX, solid_capstyle="round", zorder=3)
        ax.plot([v], [y], "o", ms=26 * PX, color=c[colour], mew=0, zorder=4)
        ax.text(hi, y, f"   {num(v, True)}  ({num(lo, True)} to {num(hi, True)})", va="center", fontsize=VALUE,
                color=c["dark"])
    if pairs:
        ax.set_yticks([k * step + 0.5 for k in range(len(pairs))], pairs)
    else:
        ax.set_yticks(ys, [r[0] for r in rows])
    style(ax, c)
    return ys


def gain_axes(ax, c: dict, lo: float, hi: float, label: str, ys: list[float]) -> None:
    """The reading rule as a band from -5 to +5, labelled once just above it, and a line at zero."""
    ax.axvspan(-5, 5, color=c["band"], lw=0, zorder=0)
    ax.axvline(0, color=c["dark"], lw=2 * PX * 2, zorder=1)
    ax.annotate("too small to count", (0, 1), xycoords=ax.get_xaxis_transform(), xytext=(0, 6 * PX),
                textcoords="offset points", ha="center", va="bottom", fontsize=TICK, color=c["dark"])
    ax.set_ylim(ys[-1] + 0.5, -0.6)
    ax.set_xlim(lo, hi)
    ax.set_xticks([x for x in range(-10, 45, 5) if lo <= x <= hi])
    ax.set_xticklabels([num(x, True).replace(".0", "") for x in range(-10, 45, 5) if lo <= x <= hi])
    ax.set_xlabel(label, fontsize=TICK, color=c["grey"], labelpad=8)


def top1_scores(c: dict, out: Path) -> None:
    s = load("replication_main.json")["strata"]["all"]["scores"]
    rows = [("Combined-score baseline", s["Lookups combined (locked)"], "strong"),
            ("Rule baseline", s["Rule-based"], "strong"),
            ("AVI", s["AVI"], "mid"), ("CADD", s["CADD"], "rest"),
            ("GPN-Star-V", s["GPN-Star-V"], "mid"), ("Evo 2", s["Evo 2"], "mid"),
            ("SpliceAI", s["SpliceAI"], "rest"), ("GPN-MSA", s["GPN-MSA"], "mid"),
            ("REVEL", s["REVEL"], "rest"), ("AlphaMissense", s["AlphaMissense"], "mid")]
    fig, ax = figure(1664, 600)
    bars(ax, rows, c, "Spiked-in variant ranked first, 500 new cases (%); thin line: 95% interval")
    save(fig, out, "v2_top1_scores", c,
         [("strong", "Baseline"), ("mid", "Foundation-model score"), ("rest", "Existing score alone")])


def indel_top1(c: dict, out: Path) -> None:
    s = load("indel_main.json")["strata"]["all"]["scores"]
    rows = [("Combined-score baseline + Evo 2", s["Lookups + Evo 2"], "mid"),
            ("Combined-score baseline", s["Lookups combined"], "strong"),
            ("Evo 2 alone", s["Evo 2"], "mid"), ("Rule baseline", s["Rule-based"], "strong"),
            ("SpliceAI", s["SpliceAI"], "rest")]
    fig, ax = figure(1728, 560)
    bars(ax, rows, c, "Spiked-in indel ranked first, 120 indel cases (%); thin line: 95% interval")
    save(fig, out, "v2_indel_top1", c,
         [("strong", "Baseline"), ("mid", "With Evo 2, alone or added"), ("rest", "Existing score alone")])


def by_group(c: dict, out: Path) -> None:
    st = load("replication_main.json")["strata"]
    fig, axes = figure(1664, 560, ncols=len(GROUPS))
    labels = ["Combined-score baseline", "Rule baseline", "Its existing score alone"]
    for ax, (key, title, tool) in zip(axes, GROUPS):
        s = st[key]["scores"]
        rows = [(pct(s["Lookups combined (locked)"]), "strong"), (pct(s["Rule-based"]), "mid"),
                (pct(s[tool]), "rest")]
        for i, ((v, lo, hi), colour) in enumerate(rows):
            ax.barh(i, v, height=0.66, color=c[colour], linewidth=0, zorder=2)
            ax.plot([lo, hi], [i, i], color=c["grey"], lw=3 * PX, solid_capstyle="round", zorder=3)
            ax.text(hi + 3, i, f"{v:.0f}", va="center", fontsize=VALUE, color=c["dark"], fontweight="bold")
        ax.set_title(f"{title} ({st[key]['cases']})\nits score: {tool}", fontsize=LABEL * 0.94, color=c["dark"],
                     pad=10, loc="left")
        ax.set_xlim(0, 135)
        ax.set_xticks([0, 50, 100])
        ax.set_ylim(len(rows) - 0.5, -0.55)
        ax.set_yticks(range(len(rows)), labels)
        style(ax, c)
    axes[0].set_xlabel("Ranked first (%)", fontsize=TICK, color=c["grey"])
    save(fig, out, "v2_by_group", c, [("strong", "Combined-score baseline"), ("mid", "Rule baseline"),
                                      ("rest", "The group's existing score alone")])


def gains_avi(c: dict, out: Path) -> None:
    rep, snv = load("replication_main.json")["strata"]["all"]["contrasts"], load("snv_main.json")["strata"]["all"]["contrasts"]
    rows = [("", pct(rep[PRIMARY]), "strong"), ("", pct(rep[SECONDARY]), "mid"),
            ("", pct(snv["Lookups + AVI minus Lookups combined"]), "strong"), ("", pct(snv[SECONDARY]), "mid")]
    fig, ax = figure(1664, 547)
    ys = dots(ax, rows, c, pairs=["500 new cases", "First 200 cases\n(cross-validated)"])
    gain_axes(ax, c, -8, 16, "Gain at top 1 from adding AVI (percentage points)", ys)
    save(fig, out, "v2_gains_avi", c, [("strong", TO_COMBINED), ("mid", TO_RULE)], marker="●")


def gains_by_group(c: dict, out: Path) -> None:
    st = load("replication_main.json")["strata"]
    rows, pairs = [], []
    for key, title in [("noncoding", "Non-coding"), ("noncoding: intronic", "Intronic"),
                       ("noncoding: splice_region", "Splice region"), ("missense", "Missense")]:
        pairs.append(f"{title}\n{st[key]['cases']} cases")
        rows.append(("", pct(st[key]["contrasts"][PRIMARY]), "strong"))
        rows.append(("", pct(st[key]["contrasts"][SECONDARY]), "mid"))
    fig, ax = figure(1664, 580)
    ys = dots(ax, rows, c, pairs=pairs)
    gain_axes(ax, c, -10, 42, "Gain at top 1 from adding AVI, 500 new cases (percentage points)", ys)
    save(fig, out, "v2_gains_by_group", c, [("strong", TO_COMBINED), ("mid", TO_RULE)], marker="●")


def gains_added(c: dict, out: Path) -> None:
    con = load("replication_main.json")["strata"]["all"]["contrasts"]
    rows = [(n, pct(con[f"Lookups + {n} (locked) minus Lookups combined (locked)"], ci="top_1_ci98.33"), "strong")
            for n in ("GPN-MSA", "GPN-Star-V", "GPN-Star-M")]
    rows.append(("Evo 2", pct(con["Lookups + Evo 2 (locked) minus Lookups combined (locked)"], ci="top_1_ci98.75"),
                 "strong"))
    rows.append(("AVI (95% interval)", pct(con[PRIMARY]), "strong"))
    fig, ax = figure(1664, 523)
    ys = dots(ax, rows, c)
    gain_axes(ax, c, -8, 12, "Gain at top 1, 500 new cases (percentage points); intervals 98.33% for GPN, "
              "98.75% for Evo 2", ys)
    save(fig, out, "v2_gains_added", c, [("strong", TO_COMBINED)], marker="●")


def first_vs_second(c: dict, out: Path) -> None:
    f200 = load("results.json", FIRST)["strata"]["all"]["added_to_lookups"]["Lookups + AlphaGenome AVI"]
    key = "Lookups + AlphaGenome AVI (locked)"
    f500 = load("results_replication.json", FIRST)["strata"]["all"]["added_to_locked_lookups"][key]
    fgrp = load("results_replication_gnomad0.001_grpmax.json", FIRST)["strata"]["all"]["added_to_locked_lookups"][key]
    v2 = load("replication_main.json")["strata"]["all"]["contrasts"][PRIMARY]
    rows = [("First run, first 200 cases", pct(f200), "rest"),
            ("First run, 500 new cases", pct(f500), "rest"),
            ("First run, 500 new, strict filter", pct(fgrp), "rest"),
            ("Second run, 500 new cases", pct(v2), "strong")]
    fig, ax = figure(1664, 523)
    ys = dots(ax, rows, c)
    gain_axes(ax, c, -8, 24, "Gain at top 1 from adding AVI to the combined-score baseline (percentage points)", ys)
    save(fig, out, "v2_first_vs_second", c, [("strong", "Second run"), ("rest", "First run")], marker="●")


def checks(c: dict, out: Path) -> None:
    rows, pairs = [], []
    for name, label in [("main", "Main filters"), ("bg1pct", "Looser filter (below 1%)"),
                        ("nomoi", "No inheritance filter"), ("all1kg", "All of 1000 Genomes left out")]:
        d = load(f"replication_{name}.json")
        con = d["strata"]["all"]["contrasts"]
        pairs.append(f"{label}\nmedian list {d['median_variants_per_case']:g}")
        rows.append(("", pct(con[PRIMARY]), "strong"))
        rows.append(("", pct(con[SECONDARY]), "mid"))
    fig, ax = figure(1664, 500)
    ys = dots(ax, rows, c, pairs=pairs)
    gain_axes(ax, c, -10, 24, "Gain at top 1 from adding AVI, 500 new cases (percentage points)", ys)
    save(fig, out, "v2_checks", c, [("strong", TO_COMBINED), ("mid", TO_RULE)], marker="●")


# Frequency sweep (PROTOCOL.md section 13, exploratory)

SWEEP_SCORES = ["AVI", "CADD", "SpliceAI", "REVEL", "AlphaMissense", "GPN-MSA", "GPN-Star-V", "GPN-Star-M", "Evo 2"]
SWEEP_BASELINES = [("Combined-score baseline", "Lookups combined (locked)"), ("Rule baseline", "Rule-based")]
SWEEP_MARKS = [(0.001, "dominant"), (0.01, "recessive"), (0.05, "BA1")]
SWEEP_GAINS = [(PRIMARY, "strong"), (SECONDARY, "mid")]


def threshold_text(af: float) -> str:
    return f"{af * 100:g}%"


def spread(ys: list[float], gap: float) -> list[float]:
    """Label heights at least gap apart, in the order of ys, moved as little as possible."""
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    out = [0.0] * len(ys)
    last = -1e9
    for i in order:
        out[i] = max(ys[i], last + gap)
        last = out[i]
    return out


def af_sweep(c: dict, out: Path) -> None:
    rows = load("replication_sweep.json")["rows"]
    afs = [r["max_af"] for r in rows.values()]
    fig, ax = figure(1664, 760)
    ends = []  # (y at the last threshold, label, colour, bold)
    for name, key in SWEEP_BASELINES:
        xs = [r["max_af"] for r in rows.values() if key in r["strata"]["all"]["scores"]]
        ys = [r["strata"]["all"]["scores"][key]["top_1"] * 100 for r in rows.values()
              if key in r["strata"]["all"]["scores"]]
        ax.plot(xs, ys, color=c["dark"], lw=3 * PX, zorder=2)
        ends.append((xs[-1], ys[-1], f"{name}  {num(ys[-1], False)}", c["dark"], True))
    for score in reversed(SWEEP_SCORES):  # AVI drawn last, on top
        pts = [(r["max_af"], *pct(r["strata"]["all"]["scores"][score])) for r in rows.values()
               if score in r["strata"]["all"]["scores"]]
        if not pts:
            continue
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        colour = c["strong"] if score == "AVI" else c["rest"]
        top = score == "AVI"
        ax.plot(xs, ys, color=colour, lw=(5 if top else 3) * PX, zorder=4 if top else 3)
        ax.plot(xs, ys, "o", ms=(22 if top else 14) * PX, color=colour, mew=0, zorder=4 if top else 3)
        if top:
            for x, v, lo, hi in pts:
                ax.plot([x, x], [lo, hi], color=c["grey"], lw=3 * PX, solid_capstyle="round", zorder=3)
                if x < afs[-1]:  # the last value is in the end label
                    ax.text(x, hi + 1.5, num(v, False), ha="center", va="bottom", fontsize=VALUE, color=c["dark"],
                            fontweight="bold", zorder=5,
                            bbox={"facecolor": c["surface"], "edgecolor": "none", "pad": 1})
        label = score + (f"  {num(ys[-1], False)}" if top else "") + (" (to 1%)" if xs[-1] < afs[-1] else "")
        ends.append((xs[-1], ys[-1], label, c["strong"] if top else c["grey"], top))
    heights = spread([e[1] for e in ends], 4.8)
    for (x, y, label, colour, bold), h in zip(ends, heights):
        ax.annotate(label, (x, y), xytext=(afs[-1] * 1.35, h), textcoords="data", va="center", ha="left",
                    fontsize=TICK, color=colour, fontweight="bold" if bold else "normal",
                    arrowprops={"arrowstyle": "-", "color": c["light"], "lw": 2 * PX, "relpos": (0, 0.5), "shrinkA": 4}
                    if abs(h - y) > 1 or x < afs[-1]
                    else None)
    for x, label in SWEEP_MARKS:
        ax.axvline(x, color=c["light"], lw=2 * PX, zorder=1)
        ax.annotate(label, (x, 1), xycoords=ax.get_xaxis_transform(), xytext=(0, 6 * PX), textcoords="offset points",
                    ha="center", va="bottom", fontsize=TICK, color=c["dark"])
    ax.set_xscale("log")
    ax.set_xlim(afs[0] / 1.6, afs[-1] * 6)
    ax.set_xticks(afs, [threshold_text(a) for a in afs])
    ax.minorticks_off()
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    style(ax, c)
    ax.tick_params(axis="y", labelsize=TICK, labelcolor=c["grey"])
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Background kept if its highest group frequency in gnomAD is below (log scale)", fontsize=TICK,
                  color=c["grey"], labelpad=8)
    ax.set_ylabel("Spiked-in variant ranked first,\n500 new cases (%)", fontsize=TICK, color=c["grey"])
    save(fig, out, "af_sweep", c, [("strong", "AVI (thin line: 95% interval)"), ("rest", "Other single scores"),
                                   ("dark", "Baselines")], marker="●", svg=True)


def af_sweep_gain(c: dict, out: Path) -> None:
    rows = load("replication_sweep.json")["rows"]
    marks = dict(SWEEP_MARKS)
    data, pairs = [], []
    for r in rows.values():
        con = r["strata"]["all"]["contrasts"]
        af = r["max_af"]
        pairs.append(f"Below {threshold_text(af)}" + (f" ({marks[af]})" if af in marks else "")
                     + f"\n{r['cases']} cases, median list {r['median_variants_per_case']:g}")
        data += [("", pct(con[key]), colour) for key, colour in SWEEP_GAINS]
    fig, ax = figure(1664, 960)
    ys = dots(ax, data, c, pairs=pairs)
    lo = min(-8, min(v[1][1] for v in data) - 2)
    hi = max(16, max(v[1][2] for v in data) + 14)
    gain_axes(ax, c, 5 * (lo // 5), 5 * -(-hi // 5), "Gain at top 1 from adding AVI, 500 new cases (percentage points)",
              ys)
    save(fig, out, "af_sweep_gain", c, [("strong", TO_COMBINED), ("mid", TO_RULE)], marker="●", svg=True)


CHARTS = (top1_scores, indel_top1, by_group, gains_avi, gains_by_group, gains_added, first_vs_second, checks)
SWEEP_CHARTS = (af_sweep, af_sweep_gain)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--brand", default="rewire", choices=BRANDS)
    p.add_argument("--out", type=Path, default=HERE / "figures" / "out")
    p.add_argument("--only", nargs="+", metavar="CHART", help="draw these charts only (function names)")
    a = p.parse_args()
    c = BRANDS[a.brand]
    a.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 22, "axes.labelcolor": c["dark"],
                         "figure.facecolor": c["surface"], "axes.facecolor": c["surface"]})
    charts = CHARTS + (SWEEP_CHARTS if (RESULTS / "replication_sweep.json").exists() else ())
    if a.only:
        charts = [ch for ch in CHARTS + SWEEP_CHARTS if ch.__name__ in a.only]
    for chart in charts:
        chart(c, a.out)
    print(f"wrote {len(charts)} charts to {a.out} ({a.brand})")


if __name__ == "__main__":
    main()
