"""Case-level ranking benchmark on the BRCA1 saturation genome editing variants.

Question: in a case with one loss-of-function variant and K other candidate variants, how often does each
score put the loss-of-function variant first, and does adding Evo 2 to existing predictors help?

Each simulated case takes one loss-of-function (LOF) variant as the "diagnostic" variant and K variants drawn
from the functional or intermediate variants at other positions. A competitor is a candidate that scores equal
to or higher than the diagnostic variant (ties count against the score). A missing score ranks lowest: a
candidate without a score never outranks, and a diagnostic variant without a score is outranked by every
candidate. Combined scores are logistic regressions on the predictors' raw scores, fitted and applied with
5-fold cross-validation grouped by genomic position, so no variant is scored by a model that saw its label or
any other allele at its position. This tests whether Evo 2 adds information to CADD and SpliceAI.

Uncertainty: 95% intervals from resampling diagnostic variants grouped by genomic position.

Inputs come from rewire-bio/evo2-cloud-operations (Evo 2 7B H100 scores, assay labels, and the CADD v1.7,
SpliceAI and REVEL extracts fetched by followup/current_predictors.py --fetch). Per-variant scores are not
written out; their licences do not allow redistribution.

  uv run --with pandas --with numpy python brca1-case-ranking/case_ranking.py --study ../evo2-cloud-operations
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261016
REPLICATES = 10      # cases per diagnostic variant
BOOTSTRAP = 1000
CASE_SIZES = [20, 100]


def load(study: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = importlib.util.spec_from_file_location("cp", study / "followup/current_predictors.py")
    cp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cp)
    variants = pd.read_csv(study / "data/derived/brca1_variants.tsv", sep="\t").set_index("id")
    raw = cp.load_scores(variants)
    scores = pd.DataFrame({"CADD v1.7": raw["cadd_v1.7"], "SpliceAI": raw["spliceai"],
                           "REVEL": raw["revel"], "Evo 2 7B": raw["evo2_7b"]})
    labels, groups = variants["lof"].to_numpy(), variants["pos"].to_numpy()
    scores["CADD + SpliceAI"] = cross_validated(scores[["CADD v1.7", "SpliceAI"]], labels, groups)
    scores["CADD + Evo 2"] = cross_validated(scores[["CADD v1.7", "Evo 2 7B"]], labels, groups)
    scores["CADD + SpliceAI + Evo 2"] = cross_validated(scores[["CADD v1.7", "SpliceAI", "Evo 2 7B"]], labels, groups)
    pct = scores[["CADD v1.7", "SpliceAI", "Evo 2 7B"]].rank(pct=True)
    scores["Average percentile CADD + SpliceAI"] = pct[["CADD v1.7", "SpliceAI"]].mean(axis=1)
    scores["Average percentile CADD + SpliceAI + Evo 2"] = pct.mean(axis=1)
    coding = (variants["region"] == "coding").astype(float).to_numpy()
    scores["Type-aware CADD + SpliceAI"] = cross_validated(type_aware(scores[["CADD v1.7", "SpliceAI"]], coding), labels, groups)
    scores["Type-aware CADD + SpliceAI + Evo 2"] = cross_validated(
        type_aware(scores[["CADD v1.7", "SpliceAI", "Evo 2 7B"]], coding), labels, groups)
    return variants, scores


def type_aware(features: pd.DataFrame, coding: np.ndarray) -> pd.DataFrame:
    """Let each score carry a different weight for coding and noncoding variants (region from annotation)."""
    out = features.copy()
    for name in features:
        out[f"{name} x coding"] = features[name].to_numpy() * coding
    out["coding"] = coding
    return out


def cross_validated(features: pd.DataFrame, labels: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Out-of-fold probabilities from a standardised logistic regression, folds grouped by position."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    x = features.to_numpy()
    assert not np.isnan(x).any(), "combined scores need complete features"
    out = np.full(len(x), np.nan)
    for train, test in GroupKFold(n_splits=5).split(x, labels, groups):
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
        model.fit(x[train], labels[train])
        out[test] = model.predict_proba(x[test])[:, 1]
    return out


def simulate(variants: pd.DataFrame, k: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Return (diagnostic index, candidate index matrix) for every simulated case."""
    lof = np.flatnonzero(variants["lof"].to_numpy() == 1)
    other = np.flatnonzero(variants["lof"].to_numpy() == 0)
    pos = variants["pos"].to_numpy()
    diag, cands = [], []
    for d in lof:
        pool = other[pos[other] != pos[d]]
        for _ in range(REPLICATES):
            diag.append(d)
            cands.append(rng.choice(pool, size=k, replace=False))
    return np.array(diag), np.array(cands)


def competitors(score: np.ndarray, diag: np.ndarray, cands: np.ndarray) -> np.ndarray:
    d = score[diag][:, None]
    c = score[cands]
    beaten = np.where(np.isnan(c), False, np.where(np.isnan(d), True, c >= d))
    return beaten.sum(axis=1)


def summarise(comp: np.ndarray, diag: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> dict:
    first = (comp == 0).astype(float)
    top5 = (comp <= 4).astype(float)
    unique = np.unique(groups[diag])
    by_group = {g: np.flatnonzero(groups[diag] == g) for g in unique}
    boots = {"first": [], "top5": []}
    for _ in range(BOOTSTRAP):
        idx = np.concatenate([by_group[g] for g in rng.choice(unique, size=len(unique), replace=True)])
        boots["first"].append(first[idx].mean())
        boots["top5"].append(top5[idx].mean())
    return {
        "ranked_first": round(first.mean(), 4), "ranked_first_ci95": [round(float(x), 4) for x in np.percentile(boots["first"], [2.5, 97.5])],
        "in_top_5": round(top5.mean(), 4), "in_top_5_ci95": [round(float(x), 4) for x in np.percentile(boots["top5"], [2.5, 97.5])],
        "median_competitors": float(np.median(comp)),
    }


PAIRS = [("CADD v1.7", "CADD + Evo 2"), ("CADD + SpliceAI", "CADD + SpliceAI + Evo 2"),
         ("CADD v1.7", "CADD + SpliceAI + Evo 2"),
         ("Type-aware CADD + SpliceAI", "Type-aware CADD + SpliceAI + Evo 2"),
         ("Average percentile CADD + SpliceAI", "Average percentile CADD + SpliceAI + Evo 2")]


def paired_difference(a: np.ndarray, b: np.ndarray, diag: np.ndarray, groups: np.ndarray, rng) -> dict:
    """b minus a for ranked first and in top 5, with position-grouped bootstrap intervals."""
    diffs = {"ranked_first": (b == 0).astype(float) - (a == 0).astype(float),
             "in_top_5": (b <= 4).astype(float) - (a <= 4).astype(float)}
    unique = np.unique(groups[diag])
    by_group = {g: np.flatnonzero(groups[diag] == g) for g in unique}
    draws = [np.concatenate([by_group[g] for g in rng.choice(unique, size=len(unique))]) for _ in range(BOOTSTRAP)]
    return {name: {"difference": round(float(d.mean()), 4),
                   "ci95": [round(float(x), 4) for x in np.percentile([d[i].mean() for i in draws], [2.5, 97.5])]}
            for name, d in diffs.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", type=Path, required=True, help="path to evo2-cloud-operations")
    parser.add_argument("--out", type=Path, default=Path(__file__).with_name("results.json"))
    args = parser.parse_args()
    variants, scores = load(args.study.resolve())
    groups = variants["pos"].to_numpy()
    coding = (variants["region"] == "coding").to_numpy()
    result = {"seed": SEED, "replicates_per_diagnostic_variant": REPLICATES, "bootstrap": BOOTSTRAP,
              "variants": len(variants), "lof_variants": int(variants["lof"].sum()),
              "coverage": {name: int(scores[name].notna().sum()) for name in scores}, "cases": {}}
    for k in CASE_SIZES:
        rng = np.random.default_rng(SEED + k)
        diag, cands = simulate(variants, k, rng)
        comp = {name: competitors(scores[name].to_numpy(), diag, cands) for name in scores}
        strata = {"all": np.ones(len(diag), bool), "coding": coding[diag], "noncoding": ~coding[diag]}
        out = {}
        for stratum, m in strata.items():
            out[stratum] = {"cases": int(m.sum()),
                            "scores": {name: summarise(c[m], diag[m], groups, rng) for name, c in comp.items()},
                            "paired": {f"{b} vs {a}": paired_difference(comp[a][m], comp[b][m], diag[m], groups, rng)
                                       for a, b in PAIRS}}
        result["cases"][f"k{k}"] = out
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    for k in CASE_SIZES:
        print(f"\nK = {k} candidates per case")
        for stratum, s in result["cases"][f"k{k}"].items():
            print(f"  {stratum} ({s['cases']} cases)")
            for name, r in s["scores"].items():
                print(f"    {name:26s} first {r['ranked_first']:.3f} {r['ranked_first_ci95']}  top5 {r['in_top_5']:.3f}  median competitors {r['median_competitors']}")
            for pair, d in s["paired"].items():
                print(f"    {pair}: first {d['ranked_first']}  top5 {d['in_top_5']}")


if __name__ == "__main__":
    main()
