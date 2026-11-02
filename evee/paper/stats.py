"""Numbers for the EVEE paper: checks on the two charts and a few values derived from EVEE's published tables.

  uv run --python 3.11 --with pandas==3.1.0 --with pyarrow==26.0.0 --with numpy==2.4.6 python evee/paper/stats.py \
      --evee ../evee-manuscript/artifacts

Run evee/evee_charts.py first: section 4 reads the SVGs it writes to evee/out/. Writes evee/paper/stats.json.

1. Pathogenic share by consequence in EVEE's full SNV set and deconfounded set, with Wilson 95% intervals computed
   here from the published counts.
2. Consequence alone as a score: each variant scored by the pathogenic share of its consequence class in that set.
   AUROC computed exactly from the counts (ties count one half), with a 95% interval from 2,000 bootstrap draws that
   resample pathogenic and benign variants separately (multinomial over classes), refitting the shares each draw.
   Also: a two-fold cross-fitted AUROC (each class's counts split at random into halves, shares fitted on one half
   and scored on the other, both ways, averaged over 200 splits); the share of pathogenic-benign pairs that fall in
   the same class; and the accuracy of predicting each class's majority label.
3. The AUROCs on the results chart (deconfounded set), with approximate 95% intervals computed here by the Hanley
   and McNeil (1982) formula from the number of pathogenic and benign variants each method scored. Where a method
   scored only part of a class, the split of the scored variants is not published; the interval then assumes the
   class's pathogenic share and is marked as assumed.
4. Check that the charts show the published values: the confounding chart's text labels and the results chart's
   marker positions are read back from the SVGs and compared with the tables.
5. The same five methods on EVEE's full SNV set, with the 95% intervals published in snv_benchmark.feather.
"""
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT_SVG = HERE.parent / "out"
SERIES = [("Covpool probe", "EVEE probe"), ("Evo2 loss", "Evo 2 zero-shot"), ("AlphaMissense", "AlphaMissense"),
          ("REVEL", "REVEL"), ("CADD v1.7", "CADD v1.7")]
ROWS = ["Missense", "Synonymous", "Intronic", "Splice Donor", "Splice Acceptor", "Nonsense"]
COLOURS = {"#b4530f": "Covpool probe", "#e8a36b": "Evo2 loss", "#2a6cb0": "AlphaMissense", "#6b4fa0": "REVEL",
           "#2e8b7a": "CADD v1.7"}
Z = 1.959964


def wilson(k, n):
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return c - h, c + h


def auc_from_counts(pos, neg, share=None):
    """AUROC of a score that gives each class its pathogenic share; pos, neg are counts per class. A share fitted
    elsewhere can be passed in."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    share = pos / np.maximum(pos + neg, 1) if share is None else np.asarray(share, float)
    order = np.argsort(share)
    pos, neg, share = pos[order], neg[order], share[order]
    num = 0.0
    for i in range(len(share)):
        below = neg[share < share[i]].sum()
        tied = neg[share == share[i]].sum()
        num += pos[i] * (below + 0.5 * tied)
    return num / (pos.sum() * neg.sum())


def consequence_only(d, rng, draws=2000):
    pos = d.n_pathogenic.to_numpy()
    neg = (d["count"] - d.n_pathogenic).to_numpy()
    est = auc_from_counts(pos, neg)
    boot = []
    for _ in range(draws):
        bp = rng.multinomial(pos.sum(), pos / pos.sum())
        bn = rng.multinomial(neg.sum(), neg / neg.sum())
        boot.append(auc_from_counts(bp, bn))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"auroc": est, "lo": lo, "hi": hi, "n": int(pos.sum() + neg.sum()), "n_pathogenic": int(pos.sum())}


def cross_fitted(d, rng, splits=200):
    """Two-fold cross-fitted AUROC of consequence alone: split each class's pathogenic and benign counts at random
    into halves, fit the shares on one half and score the other, both ways; mean over the splits."""
    pos = d.n_pathogenic.to_numpy()
    neg = (d["count"] - d.n_pathogenic).to_numpy()
    vals = []
    for _ in range(splits):
        pa, na = rng.binomial(pos, 0.5), rng.binomial(neg, 0.5)
        pb, nb = pos - pa, neg - na
        sa, sb = pa / np.maximum(pa + na, 1), pb / np.maximum(pb + nb, 1)
        vals.append((auc_from_counts(pb, nb, sa) + auc_from_counts(pa, na, sb)) / 2)
    return {"auroc": float(np.mean(vals)), "splits": splits}


def pairs_and_accuracy(d):
    """Share of pathogenic-benign pairs in the same class, and accuracy of predicting each class's majority label."""
    pos = d.n_pathogenic.to_numpy().astype(float)
    neg = (d["count"] - d.n_pathogenic).to_numpy().astype(float)
    return {"within_class_pair_share": float((pos * neg).sum() / (pos.sum() * neg.sum())),
            "majority_rule_accuracy": float(np.maximum(pos, neg).sum() / (pos + neg).sum())}


def hanley_mcneil(a, n1, n0):
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    se = math.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))
    return max(0.0, a - Z * se), min(1.0, a + Z * se)


def svg_results(path):
    """Read the results chart's markers back to AUROC values, by colour, top to bottom."""
    s = path.read_text()
    uses = re.findall(r'<use xlink:href="#(\w+)" x="([\d.]+)" y="([\d.]+)" style="fill: (#\w+); stroke: #fdfcfa', s)
    ticks = sorted({float(x) for _, x, _, _ in
                    re.findall(r'<use xlink:href="#(\w+)" x="([\d.]+)" y="([\d.]+)" style="fill: (#5e6874)', s)})
    x0, x1 = ticks[0], ticks[-1]  # 0.70 and 1.00
    got = {}
    for _, x, y, fill in uses:
        got.setdefault(COLOURS[fill], []).append((float(y), 0.7 + (float(x) - x0) / (x1 - x0) * 0.3))
    return {m: [v for _, v in sorted(pts)] for m, pts in got.items()}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--evee", type=Path, required=True)
    a = p.parse_args()
    rng = np.random.default_rng(20261102)
    rng_cf = np.random.default_rng(20261103)  # separate stream, so the bootstrap intervals do not change
    out = {}

    # 1 and 2
    dc = pd.read_feather(a.evee / "dataset_characterization.feather")
    out["composition"] = {}
    out["consequence_only"] = {}
    for key, prefix in (("full", "Naive"), ("deconfounded", "Deconf")):
        d = dc[dc.dataset.str.startswith(prefix)]
        rows = []
        for _, r in d.iterrows():
            lo, hi = wilson(r.n_pathogenic, r["count"])
            rows.append({"consequence": r.consequence, "count": int(r["count"]), "n_pathogenic": int(r.n_pathogenic),
                         "share_pct": 100 * r.n_pathogenic / r["count"], "lo_pct": 100 * lo, "hi_pct": 100 * hi})
        out["composition"][key] = {"total": int(d["count"].sum()), "pathogenic": int(d.n_pathogenic.sum()),
                                   "intronic_plus_synonymous": int(d[d.consequence.isin(["Intronic", "Synonymous"])]
                                                                   ["count"].sum()),
                                   "rows": rows}
        out["consequence_only"][key] = consequence_only(d, rng)
        out["consequence_only"][key]["cross_fitted"] = cross_fitted(d, rng_cf)
        out["consequence_only"][key].update(pairs_and_accuracy(d))

    # 3
    db = pd.read_feather(a.evee / "deconf_benchmark.feather")
    chart = []
    for method, label in SERIES:
        sub = db[db.method == method].set_index("consequence")
        for row in ["Overall"] + ROWS:
            if row not in sub.index:
                continue
            r = sub.loc[row]
            n = int(r.n_valid)
            assumed = n != int(r.n_total)
            n1 = int(r.n_pathogenic) if not assumed else round(n * r.n_pathogenic / r.n_total)
            if method == "AlphaMissense" and row == "Missense":
                # AlphaMissense scores exactly the canonical missense variants: same n, published split
                can = sub.loc["Missense (canonical)"]
                if int(can.n_total) == n:
                    n1, assumed = int(can.n_pathogenic), False
            lo, hi = hanley_mcneil(float(r.auroc), n1, n - n1)
            chart.append({"method": method, "label": label, "consequence": row, "n_scored": n,
                          "n_pathogenic_scored": n1, "split_assumed": assumed, "auroc": float(r.auroc),
                          "lo": lo, "hi": hi, "on_chart": row != "Overall" and n >= 100})
    out["deconfounded_chart"] = chart

    # 4
    check = {"results_chart": [], "confounding_chart": []}
    svg = svg_results(OUT_SVG / "evee_results.svg")
    for method, _ in SERIES:
        want = [c["auroc"] for c in chart if c["method"] == method and c["on_chart"]]
        want_rows = [c["consequence"] for c in chart if c["method"] == method and c["on_chart"]]
        got = svg.get(method, [])
        for row, w, g in zip(want_rows, want, got):
            check["results_chart"].append({"method": method, "consequence": row, "table": w, "chart": g,
                                           "diff": abs(w - g)})
        if len(got) != len(want):
            check["results_chart"].append({"method": method, "error": f"{len(got)} markers, {len(want)} expected"})
    s = (OUT_SVG / "evee_confounding.svg").read_text()
    labels = re.findall(r"([\d.]+)%   \(n = ([\d,]+)\)", s)
    full = {r["count"]: r for r in out["composition"]["full"]["rows"]}
    for pct, n in labels:
        r = full[int(n.replace(",", ""))]
        check["confounding_chart"].append({"consequence": r["consequence"], "chart_pct": float(pct),
                                           "table_pct": round(r["share_pct"], 1),
                                           "match": f"{r['share_pct']:.1f}" == pct})
    check["max_results_diff"] = max(c.get("diff", 1) for c in check["results_chart"])
    check["all_confounding_match"] = all(c["match"] for c in check["confounding_chart"])
    out["chart_check"] = check

    # 5
    sb = pd.read_feather(a.evee / "snv_benchmark.feather")
    full_rows = []
    for method, label in SERIES:
        sub = sb[sb.method == method].set_index("consequence")
        for row in ["Overall"] + ROWS:
            if row in sub.index:
                r = sub.loc[row]
                full_rows.append({"method": method, "label": label, "consequence": row, "n_scored": int(r.n_valid),
                                  "auroc": float(r.auroc), "lo": None if pd.isna(r.auroc_lo) else float(r.auroc_lo),
                                  "hi": None if pd.isna(r.auroc_hi) else float(r.auroc_hi)})
    out["full_set"] = full_rows

    (HERE / "stats.json").write_text(json.dumps(out, indent=1))
    print("wrote", HERE / "stats.json")
    print("consequence only:", {k: round(v["auroc"], 4) for k, v in out["consequence_only"].items()})
    print("chart check: max AUROC difference", round(check["max_results_diff"], 5),
          "; confounding labels all match:", check["all_confounding_match"])


if __name__ == "__main__":
    main()
