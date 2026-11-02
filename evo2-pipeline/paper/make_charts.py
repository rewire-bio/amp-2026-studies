"""The paper's three charts, from run/summary.json and run/ranks.tsv (written by record_run.sh).

  uv run -q --python 3.12 --with matplotlib==3.11.2 python make_charts.py --out figures/charts

  case_ranks     where each score ranks the test case's diagnostic variant among its 100 variants, with the range
                 of ranks a tied score could take
  h100_vs_l40s   each variant's Evo 2 rank from the H100 scores against its rank from the L40S scores
  cloud_runs     the recorded cloud runs whose scores the replay profile uses: time from launch to workload start,
                 the workloads, and shutdown, with the cost of each run

Each chart is written as SVG and PNG.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
# Same colour-blind safe slots as the other studies' charts, plus grey and ink for text.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GREY, LIGHT, INK, MUTED = "#8a8985", "#d9d8d4", "#0b0b0b", "#52514e"


def style() -> None:
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "figure.dpi": 150,
                         "font.family": ["Helvetica", "Arial", "DejaVu Sans"]})


def save(fig, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out / name}.svg and .png")


def case_ranks(summary: dict, rows: list[dict], out: Path) -> None:
    d = summary["diagnostic"]
    n = summary["case"]["variants"]
    series = [("CADD v1.7", "cadd"), ("SpliceAI", "spliceai"), ("Evo 2 7B, H100", "evo2_h100"),
              ("Evo 2 7B, L40S", "evo2_l40s")]
    fig, ax = plt.subplots(figsize=(6.4, 2.3))
    for i, (label, key) in enumerate(series):
        y = len(series) - 1 - i
        ax.plot([1, n], [y, y], color=LIGHT, lw=1, zorder=1)
        lo, hi, rank = d[key]["rank_min"], d[key]["rank_max"], d[key]["rank"]
        if hi > lo:
            ax.plot([lo, hi], [y, y], color=BLUE, lw=3, alpha=0.35, solid_capstyle="round", zorder=2)
            ax.text((lo + hi) / 2, y + 0.28, f"tied with {hi - lo} others: ranks {lo} to {hi}", ha="center",
                    va="bottom", fontsize=8, color=MUTED)
        ax.scatter([rank], [y], s=48, color=BLUE, edgecolor="white", linewidth=1.5, zorder=3)
        ax.text(rank, y - 0.3, f"{rank}", ha="center", va="top", fontsize=8, color=INK)
    ax.set_yticks(range(len(series)), [s[0] for s in reversed(series)])
    ax.set_xlim(0, n + 1)
    ax.set_ylim(-0.7, len(series) - 0.3)
    ax.set_xticks([1, 10, 20, 40, 60, 80, 100])
    ax.set_xlabel(f"Rank of the diagnostic variant among the case's {n} variants (1 = top)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    save(fig, out, "case_ranks")


def h100_vs_l40s(summary: dict, rows: list[dict], out: Path) -> None:
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    n = summary["case"]["variants"]
    ax.plot([1, n], [1, n], color=LIGHT, lw=1, zorder=1)
    others = [r for r in rows if r["diagnostic"] == "0"]
    diag = [r for r in rows if r["diagnostic"] == "1"][0]
    ax.scatter([int(r["evo2_rank_h100"]) for r in others], [int(r["evo2_rank_l40s"]) for r in others], s=14,
               color=GREY, edgecolor="white", linewidth=0.6, zorder=2, label="Other variants")
    ax.scatter([int(diag["evo2_rank_h100"])], [int(diag["evo2_rank_l40s"])], s=48, color=BLUE, edgecolor="white",
               linewidth=1.5, zorder=3, label="Diagnostic variant")
    ax.annotate(f"diagnostic: {diag['evo2_rank_h100']} on H100, {diag['evo2_rank_l40s']} on L40S",
                (int(diag["evo2_rank_h100"]), int(diag["evo2_rank_l40s"])), xytext=(40, 24), textcoords="data",
                fontsize=8, color=INK, va="center", arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
    m = summary["l40s_vs_h100"]
    ax.text(n, 4, f"{m['variants_whose_evo2_rank_changes']} of {n} ranks change\nmedian move "
            f"{m['median_rank_move']:.0f}, largest {m['largest_rank_move']}",
            ha="right", va="bottom", fontsize=8, color=MUTED)
    ax.set_xlim(0, n + 2)
    ax.set_ylim(0, n + 2)
    ax.set_xlabel("Evo 2 rank from H100 scores")
    ax.set_ylabel("Evo 2 rank from L40S scores")
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    save(fig, out, "h100_vs_l40s")


def cloud_runs(summary: dict, out: Path) -> None:
    runs = summary["cloud_run_recorded_b75c67aa"]
    labels = {"P2": "H100, Google Cloud\nspot (scores replayed)", "P3": "H100, Runpod\non demand",
              "P5": "L40S, Runpod\non demand"}
    t = lambda s: datetime.fromisoformat(s.replace("Z", "+00:00"))  # noqa: E731
    fig, ax = plt.subplots(figsize=(6.4, 2.2))
    order = ["P2", "P3", "P5"]
    for i, pid in enumerate(order):
        r = runs[pid]
        e = r["events_utc"]
        start = t(e["launch_requested"])
        segs = [((t(e["workload_start"]) - start).total_seconds() / 60, GREY, "Launch to workload start"),
                ((t(e["workload_end"]) - t(e["workload_start"])).total_seconds() / 60, BLUE, "Workloads"),
                ((t(e["terminated"]) - t(e["workload_end"])).total_seconds() / 60, LIGHT, "Shutdown")]
        y = len(order) - 1 - i
        left = 0.0
        for width, colour, name in segs:
            ax.barh(y, width, left=left, color=colour, height=0.55, edgecolor="white", linewidth=1,
                    label=name if i == 0 else None)
            left += width
        ax.text(left + 2, y, f"{left:.0f} min, USD {r['cost']['usd']:.2f}", va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(order)), [labels[p] for p in reversed(order)])
    ax.set_xlabel("Minutes from launch request")
    ax.set_xlim(0, 200)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.45, 1.0), fontsize=8)
    save(fig, out, "cloud_runs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "figures/charts")
    args = parser.parse_args()
    style()
    summary = json.loads((HERE / "run/summary.json").read_text())
    with (HERE / "run/ranks.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    case_ranks(summary, rows, args.out)
    h100_vs_l40s(summary, rows, args.out)
    cloud_runs(summary, args.out)


if __name__ == "__main__":
    main()
