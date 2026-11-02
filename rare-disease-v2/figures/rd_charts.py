"""Dot-and-interval charts of the second run's results (PROTOCOL.md section 11), from case_benchmark.py's JSON.

Four charts, each written as SVG and PNG:

  gains_primary      the primary contrast (Lookups + AVI minus Lookups combined) and the key secondary (Rule + AVI
                     minus Rule-based), top 1, for each study with results; dashed lines at -5, 0 and +5 points
  top1_by_group      top 1 by diagnostic group, and by class within the non-coding group, for the single tools,
                     Lookups combined, Rule-based, Rule + AVI and Lookups + AVI; panels under 20 cases are marked small
  gains_added        the top-1 gain from every added score (Lookups + X minus Lookups combined, and Rule + AVI minus
                     Rule-based); the GPN family on the 500 with its Bonferroni 98.33% interval
  gains_gnomad       the primary and key secondary in cases whose diagnostic is found in full gnomAD and in cases
                     whose diagnostic is absent, and the difference absent minus found (section 9, shared draws)

A study with results but no AVI contrast (the indel study: AVI scores SNVs only) gets a row saying so.

Gains are in percentage points. Replication rows use the locked models; snv and indel rows are cross-validated and
labelled so. Reads common.RESULTS/<study>_<analysis>.json for snv, indel and replication.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python figures/rd_charts.py \\
      [--analysis main] [--out figures/out]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import common  # noqa: E402

# Colour-blind safe slots (validated categorical palette, first three slots), plus grey and ink for text.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GREY, INK, MUTED = "#8a8985", "#0b0b0b", "#52514e"
STUDIES = [("replication", "500 replication cases, locked"), ("snv", "200 SNV cases, cross-validated"),
           ("indel", "120 indel cases, cross-validated")]
LOCKED = " (locked)"
LOOKUPS, PLUS_AVI, RULE, RULE_AVI = "Lookups combined", "Lookups + AVI", "Rule-based", "Rule + AVI"
SINGLES = ["CADD", "SpliceAI", "REVEL", "AlphaMissense", "AVI", "GPN-MSA", "GPN-Star-V", "GPN-Star-M", "Evo 2"]


def style() -> None:
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "figure.dpi": 150})


def load(analysis: str) -> dict[str, dict]:
    out = {}
    for study, _ in STUDIES:
        f = common.RESULTS / f"{study}_{analysis}.json"
        if f.exists():
            out[study] = json.loads(f.read_text())
    if not out:
        raise SystemExit(f"no {analysis} results in {common.RESULTS}")
    return out


def fit(study: str) -> str:
    return LOCKED if study == "replication" else ""


def save(fig, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out / name}.svg and .png")


def dot_rows(ax, rows: list[tuple], xlabel: str, guides: tuple = ()) -> None:
    """rows: (label, estimate, [low, high], colour). Top row first. Values in percentage points.

    A row whose estimate is None is a note: its colour slot holds the text, written across the panel.
    """
    for x in guides:
        ax.axvline(x, color=MUTED if x == 0 else GREY, lw=0.8, ls="--", zorder=0)
    for i, (label, est, ci, colour) in enumerate(rows):
        y = len(rows) - 1 - i
        if est is None:
            ax.text(0.02, y, colour, transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=MUTED,
                    style="italic", zorder=4, bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})
            continue
        ax.plot(ci, [y, y], color=colour, lw=2, solid_capstyle="round")
        ax.plot(est, y, "o", ms=7, color=colour, mec="white", mew=1.5, zorder=3)
        ax.annotate(f"{est:+.1f}" if guides else f"{est:.0f}", (ci[1], y), xytext=(5, 0),
                    textcoords="offset points", va="center", fontsize=8, color=MUTED)
    ax.set_yticks(range(len(rows)), [r[0] for r in reversed(rows)])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#e6e5e1", lw=0.6)
    ax.set_axisbelow(True)


def pts(row: dict, key: str = "top_1", ci: str = "top_1_ci95") -> tuple[float, list[float]]:
    return 100 * row[key], [100 * v for v in row[ci]]


def no_avi(study: str) -> str:
    why = "AVI scores SNVs only" if study == "indel" else "AVI has no coverage in these results"
    return f"No AVI contrast in this study ({why})"


def gains_primary(res: dict, out: Path) -> None:
    rows = []
    for study, label in STUDIES:
        if study not in res:
            continue
        c = res[study]["strata"]["all"].get("contrasts", {})
        f = fit(study)
        before = len(rows)
        for name, colour, short in ((f"{PLUS_AVI}{f} minus {LOOKUPS}{f}", BLUE, "Lookups + AVI vs Lookups combined"),
                                    (f"{RULE_AVI} minus {RULE}", ORANGE, "Rule + AVI vs Rule-based")):
            if name in c:
                est, ci = pts(c[name])
                rows.append((f"{short}\n{label}", est, ci, colour))
        if len(rows) == before:
            rows.append((label, None, None, no_avi(study)))
    if all(r[1] is None for r in rows):
        print("gains_primary: no primary or key secondary contrast in the results; skipped")
        return
    fig, ax = plt.subplots(figsize=(7, 0.6 * len(rows) + 1.2))
    dot_rows(ax, rows, "Top-1 gain, percentage points (95% paired interval)", guides=(-5, 0, 5))
    ax.set_title("Gain from AVI, top 1 (dashed: no gain, and plus or minus 5 points)", loc="left", fontsize=10)
    save(fig, out, "gains_primary")


def panels(strata: dict) -> list[str]:
    """All cases, then each diagnostic group followed by its class panels (strata named "group: class")."""
    groups = [g for g in strata if g != "all" and ":" not in g]
    out = ["all"]
    for g in groups:
        out += [g] + [s for s in strata if s.startswith(f"{g}: ")]
    return out


def top1_by_group(res: dict, out: Path) -> None:
    for study, label in STUDIES:
        if study not in res:
            continue
        strata = res[study]["strata"]
        f = fit(study)
        wanted = SINGLES + [RULE, RULE_AVI, f"{LOOKUPS}{f}", f"{PLUS_AVI}{f}"]
        groups = panels(strata)
        scores = [s for s in wanted if s in strata["all"]["scores"]]
        ncols = min(len(groups), 5)
        nrows = math.ceil(len(groups) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(2.3 * ncols + 1.5, (0.32 * len(scores) + 1.2) * nrows + 0.4),
                                 sharey=True, sharex=True, squeeze=False, layout="constrained")
        axes = axes.ravel()
        for i, (ax, g) in enumerate(zip(axes, groups)):
            rows = []
            for s in scores:
                colour = (ORANGE if s.startswith(PLUS_AVI) or s == RULE_AVI else
                          BLUE if s in (RULE, f"{LOOKUPS}{f}") else GREY)
                est, ci = pts(strata[g]["scores"][s])
                rows.append((s.replace(LOCKED, ", locked"), est, ci, colour))
            dot_rows(ax, rows, "Top 1, %")
            small = ", small (under 20)" if strata[g].get("small") else ""
            ax.set_title(f"{g.replace('_', ' ')}\n{strata[g]['cases']} cases{small}", fontsize=9, loc="left",
                         color=MUTED if small else INK)
            ax.set_xlim(0, 105)
            if i < len(groups) - ncols:
                ax.set_xlabel("")
        for i, ax in enumerate(axes[:len(groups)]):
            if i + ncols >= len(groups):  # no panel below: show its own x tick labels
                ax.tick_params(axis="x", labelbottom=True)
                ax.set_xlabel("Top 1, %")
        for ax in axes[len(groups):]:
            ax.set_visible(False)
        fig.suptitle(f"Top 1 by diagnostic group, {label} (95% intervals)", x=0.01, ha="left", fontsize=10)
        save(fig, out, f"top1_by_group_{study}")


def gains_added(res: dict, out: Path) -> None:
    rows = []
    for study, label in STUDIES:
        if study not in res:
            continue
        c = res[study]["strata"]["all"].get("contrasts", {})
        f = fit(study)
        for x in ["AVI", "GPN-MSA", "GPN-Star-V", "GPN-Star-M", "Evo 2"]:
            name = f"Lookups + {x}{f} minus {LOOKUPS}{f}"
            if name not in c:
                continue
            adjusted = "top_1_ci98.33" in c[name] and study == "replication"
            est, ci = pts(c[name], ci="top_1_ci98.33" if adjusted else "top_1_ci95")
            note = ", 98.33% interval" if adjusted else ""
            rows.append((f"{x} added to lookups\n{label}{note}", est, ci, BLUE if x == "AVI" else AQUA))
        name = f"{RULE_AVI} minus {RULE}"
        if name in c:
            est, ci = pts(c[name])
            rows.append((f"AVI added to the rule\n{label}", est, ci, ORANGE))
        if f"Lookups + AVI{f} minus {LOOKUPS}{f}" not in c and name not in c:
            rows.append((label, None, None, no_avi(study)))
    if all(r[1] is None for r in rows):
        print("gains_added: no added-score contrasts in the results; skipped")
        return
    fig, ax = plt.subplots(figsize=(7, 0.55 * len(rows) + 1.2))
    dot_rows(ax, rows, "Top-1 gain, percentage points (95% paired interval unless marked)", guides=(-5, 0, 5))
    ax.set_title("Top-1 gain from each added score", loc="left", fontsize=10)
    save(fig, out, "gains_added")


def gains_gnomad(res: dict, out: Path) -> None:
    """Section 9: primary and key secondary in found and absent cases, and absent minus found, on shared draws."""
    rows = []
    for study, label in STUDIES:
        if study not in res:
            continue
        split = res[study].get("gnomad_split")
        if not split:
            continue
        f = fit(study)
        before = len(rows)
        for name, short, colour in ((f"{PLUS_AVI}{f} minus {LOOKUPS}{f}", "Lookups + AVI vs Lookups combined", BLUE),
                                    (f"{RULE_AVI} minus {RULE}", "Rule + AVI vs Rule-based", ORANGE)):
            for stratum, word in (("gnomad: found", "found in"), ("gnomad: absent", "absent from")):
                c = split.get(stratum, {}).get("contrasts", {})
                if name in c:
                    est, ci = pts(c[name])
                    rows.append((f"{short}, {word} gnomAD ({split[stratum]['cases']} cases)\n{label.split(' ', 1)[1]}", est, ci,
                                 colour))
            d = split.get("difference", {}).get(name)
            if d:
                est, ci = pts(d, key="top_1_absent_minus_found")
                flag = ", possibly gnomAD presence" if d.get("possibly_gnomad_presence") else ""
                rows.append((f"{short}, absent minus found{flag}\n{label}", est, ci, GREY))
        if len(rows) == before:
            rows.append((label, None, None, no_avi(study)))
    if all(r[1] is None for r in rows):
        print("gains_gnomad: no gnomAD split contrasts in the results; skipped")
        return
    fig, ax = plt.subplots(figsize=(7.5, 0.6 * len(rows) + 1.2), layout="constrained")
    dot_rows(ax, rows, "Top-1 gain, percentage points (95% paired interval, shared draws)", guides=(-5, 0, 5))
    ax.set_title("Gain from AVI by whether the diagnostic is in full gnomAD", loc="left", fontsize=10)
    save(fig, out, "gains_gnomad")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--analysis", default="main")
    p.add_argument("--out", type=Path, default=HERE / "figures/out")
    a = p.parse_args()
    style()
    res = load(a.analysis)
    gains_primary(res, a.out)
    top1_by_group(res, a.out)
    gains_added(res, a.out)
    gains_gnomad(res, a.out)


if __name__ == "__main__":
    main()
