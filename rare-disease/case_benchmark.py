"""Rare-disease case benchmark: where does each score rank the diagnostic variant within a real background?

Each case (build_cases.py) is one 1000 Genomes individual's rare SNVs in intellectual-disability panel regions
plus one ClinVar P/LP diagnostic SNV. For every score, count the background variants that score equal to or
higher than the diagnostic variant (ties count against the score). A missing score ranks lowest, as a lookup
that returns nothing would in a pipeline that sorts by that score.

Combined scores are logistic regressions fitted with 5-fold cross-validation grouped by 1000 Genomes sample
(four cases share each background), so no case is scored by a model that saw it. Inputs: each score (missing
set to 0, with a missing indicator), the variant's consequence class from CADD, and score-by-class terms, so a
score can count for missense but not for intronic variants, as a pipeline would weight it.

Scores: established lookups (CADD v1.7 PHRED, SpliceAI max delta, REVEL, AlphaMissense); foundation-model
tables (GPN-MSA, GPN-Star vertebrate and mammal, AlphaGenome AVI when scores_avi.tsv exists); and Evo 2 7B scored
by us (minus delta, so higher means more damaging). Each foundation-model score is also added to the combined
lookups, one at a time, to measure what it adds.

Consequence class: CADD's for SNVs (CADD covers every SNV); annotate_indels.py for indels, because CADD's indel
table covers only indels seen in gnomAD.

Replication (--study replication, 500 new SNV cases in data/replication/): the combined models are fitted as above
on the new cases alone, and also locked: fitted once on all 200 SNV cases and applied unchanged to the new cases.
Evo 2 was not run on the new cases, so it is left out of both.

Frequency filter (--gnomad-max-af, from gnomad_af.py): drop every variant whose gnomAD v4.1 genome AF (or, with
--grpmax, AF_grpmax) is at or above the threshold, as a diagnostic pipeline would before ranking, and drop the cases
whose diagnostic variant is itself removed. The locked replication models are then fitted on the filtered 200 cases.

  python rare-disease/case_benchmark.py [--study snv|indel|cancer|replication]   # writes results[_<study>].json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
# Evo 2 runs: run/ (SNV study, plus scores_extra.tsv), run3/ (indel and cancer windows not scored before).
EVO2 = [HERE / "evo2/run", HERE / "evo2/run3"]
FM_TABLES = ["GPN-MSA", "GPN-Star-V", "GPN-Star-M", "AlphaGenome AVI"]
SEED = 20261016
BOOTSTRAP = 1000
KEY = ["chrom", "pos", "ref", "alt"]
TOP = [1, 5, 10]


def load(data: Path, evo2: bool = True) -> tuple[pd.DataFrame, list[str]]:
    cases = pd.read_csv(data / "cases.tsv", sep="\t")
    sources = {"scores_cadd.tsv": {"cadd_phred": "CADD v1.7", "spliceai_max": "SpliceAI"},
               "scores_revel.tsv": {"revel": "REVEL"}, "scores_am.tsv": {"alphamissense": "AlphaMissense"},
               "scores_gpnmsa.tsv": {"gpn_msa": "GPN-MSA"}, "scores_gpnstar_v.tsv": {"gpnstar_v": "GPN-Star-V"},
               "scores_gpnstar_m.tsv": {"gpnstar_m": "GPN-Star-M"}, "scores_avi.tsv": {"avi": "AlphaGenome AVI"}}
    singles = []
    for name, cols in sources.items():
        path = data / name
        if not path.exists():
            continue
        extra = ["cadd_consequence"] if name == "scores_cadd.tsv" else []
        s = pd.read_csv(path, sep="\t")[KEY + list(cols) + extra].rename(columns=cols)
        cases = cases.merge(s, on=KEY, how="left")
        singles += list(cols.values())
    files = sorted(f for d in EVO2 for f in d.glob("scores*.tsv")) if evo2 else []
    if files:
        e = pd.concat([pd.read_csv(f, sep="\t") for f in files]).drop_duplicates(KEY)[KEY + ["evo2_delta"]]
        e["Evo 2 7B"] = -e["evo2_delta"]
        cases = cases.merge(e[KEY + ["Evo 2 7B"]], on=KEY, how="left")
        singles.append("Evo 2 7B")
    cls = cases["cadd_consequence"].map(CLASSES).fillna("other")
    indel_file = data / "consequence_indels.tsv"
    if indel_file.exists():
        ind = cases[KEY].merge(pd.read_csv(indel_file, sep="\t"), on=KEY, how="left")["indel_class"]
        cls = cls.where(ind.isna().to_numpy(), ind.to_numpy())
    cases["class"] = cls
    cases["is_snv"] = (cases["ref"].str.len() == 1) & (cases["alt"].str.len() == 1)
    # Labs run SpliceAI on indels missing from its precomputed files; when that has been done (run_spliceai.sh),
    # indels take the run score, so SpliceAI covers every indel equally.
    run = data / "scores_spliceai_run.tsv"
    cases["spliceai_is_run"] = False
    if run.exists():
        r = cases[KEY].merge(pd.read_csv(run, sep="\t"), on=KEY, how="left")["spliceai_run"]
        indel = ~cases["is_snv"]
        cases.loc[indel, "SpliceAI"] = r[indel].to_numpy()
        cases.loc[indel, "spliceai_is_run"] = True
    return cases, singles


# Consequence classes from CADD's most severe annotation. A pipeline always knows the consequence, so the
# combined models may weight each score differently by class.
CLASSES = {"STOP_GAINED": "lof", "STOP_LOST": "lof", "CANONICAL_SPLICE": "splice", "SPLICE_SITE": "splice",
           "NON_SYNONYMOUS": "missense", "SYNONYMOUS": "synonymous", "5PRIME_UTR": "utr", "3PRIME_UTR": "utr",
           "INTRONIC": "intronic"}
INDEL_CLASSES = ["frameshift", "inframe"]  # splice, utr and intronic indels share the SNV classes


# Lookup tables whose indel coverage depends on gnomAD: 1000 Genomes is part of gnomAD, so background indels are
# almost always present and ClinVar diagnostic indels mostly absent. Coverage would give the answer away, so the
# combined models treat these scores as unavailable for every indel.
GNOMAD_INDEL_TABLES = ["CADD v1.7", "SpliceAI"]


def features(cases: pd.DataFrame, names: list[str]) -> pd.DataFrame:
    """Scores (missing set to 0), a missing flag on SNVs only, consequence class, and score x class interactions.

    Missing flags are SNV-only: for indels, whether a table has a score reflects the study design (every indel
    table here is incomplete or gnomAD-based), not anything a pipeline could use fairly.
    """
    out = pd.DataFrame(index=cases.index)
    cls = cases["class"]
    snv = cases["is_snv"]
    for n in names:
        unfair = pd.Series(n in GNOMAD_INDEL_TABLES, index=cases.index)
        if n == "SpliceAI":
            unfair &= ~cases["spliceai_is_run"]
        values = cases[n].where(snv | ~unfair)
        missing = values.isna() & snv
        out[n] = values.fillna(0.0)
        if missing.any():
            out[f"{n} missing"] = missing.astype(float)
    for c in sorted(set(CLASSES.values()) | set(INDEL_CLASSES)):
        if not (cls == c).any():
            continue
        is_c = (cls == c).astype(float)
        out[f"is {c}"] = is_c
        for n in names:
            out[f"{n} x {c}"] = out[n] * is_c
    return out


def gnomad_filter(cases: pd.DataFrame, max_af: float | None, grpmax: bool) -> pd.DataFrame:
    """Drop variants at or above max_af in gnomAD, then the cases that lost their diagnostic variant."""
    if max_af is None:
        return cases
    col = "gnomad_af_grpmax" if grpmax else "gnomad_af"
    af = cases[KEY].merge(pd.read_csv(DATA / "gnomad_af.tsv", sep="\t"), on=KEY, how="left")[col]
    if af.isna().any():
        raise SystemExit(f"{int(af.isna().sum())} variants have no gnomad_af.tsv row; rerun gnomad_af.py")
    kept = cases[(af < max_af).to_numpy()]
    solved = kept.loc[kept["diagnostic"] == 1, "case"]
    return kept[kept["case"].isin(solved)].reset_index(drop=True)


def classifier():
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced"))


def cross_validated(x: pd.DataFrame, labels: np.ndarray, groups: np.ndarray) -> np.ndarray:
    from sklearn.model_selection import GroupKFold

    x = x.to_numpy()
    out = np.full(len(x), np.nan)
    for train, test in GroupKFold(n_splits=5).split(x, labels, groups):
        model = classifier()
        model.fit(x[train], labels[train])
        out[test] = model.predict_proba(x[test])[:, 1]
    return out


def locked(train: pd.DataFrame, test: pd.DataFrame, names: list[str]) -> np.ndarray:
    """Fit on every training case, then score the test cases with that model unchanged.

    Feature columns follow the training cases: a missing flag or class the training cases lack is dropped.
    """
    x = features(train, names)
    model = classifier().fit(x.to_numpy(), train["diagnostic"].to_numpy())
    return model.predict_proba(features(test, names).reindex(columns=x.columns, fill_value=0.0).to_numpy())[:, 1]


def competitors(cases: pd.DataFrame, score: str) -> pd.Series:
    """Background variants scoring at or above the diagnostic variant, per case."""
    def one(g: pd.DataFrame) -> int:
        d = g.loc[g["diagnostic"] == 1, score].iloc[0]
        bg = g.loc[g["diagnostic"] == 0, score].dropna()
        return len(g) - 1 if np.isnan(d) else int((bg >= d).sum())
    return cases.groupby("case")[[score, "diagnostic"]].apply(one)


def summarise(comp: pd.Series, sample: pd.Series, rng) -> dict:
    samples = sample.unique()
    idx = {s: np.flatnonzero(sample.to_numpy() == s) for s in samples}
    draws = [np.concatenate([idx[s] for s in rng.choice(samples, size=len(samples))]) for _ in range(BOOTSTRAP)]
    c = comp.to_numpy()
    out = {"median_competitors": float(np.median(c))}
    for k in TOP:
        hit = (c < k).astype(float)
        out[f"top_{k}"] = round(hit.mean(), 4)
        out[f"top_{k}_ci95"] = [round(float(x), 4) for x in np.percentile([hit[i].mean() for i in draws], [2.5, 97.5])]
    return out


def paired(comp_a: pd.Series, comp_b: pd.Series, sample: pd.Series, rng) -> dict:
    """Difference in top-k rate, a minus b, on the same cases, with a bootstrap by sample."""
    samples = sample.unique()
    idx = {s: np.flatnonzero(sample.to_numpy() == s) for s in samples}
    draws = [np.concatenate([idx[s] for s in rng.choice(samples, size=len(samples))]) for _ in range(BOOTSTRAP)]
    a, b = comp_a.to_numpy(), comp_b.loc[comp_a.index].to_numpy()
    out = {}
    for k in TOP:
        d = (a < k).astype(float) - (b < k).astype(float)
        out[f"top_{k}"] = round(d.mean(), 4)
        out[f"top_{k}_ci95"] = [round(float(x), 4) for x in np.percentile([d[i].mean() for i in draws], [2.5, 97.5])]
    return out


# ClinGen PP3 thresholds (Bergquist et al., Genet Med 2025, Table 1): Supporting, Moderate, Strong.
PP3 = {"AlphaMissense": [0.792, 0.906, 0.990], "REVEL": [0.644, 0.773, 0.932]}


def pp3(cases: pd.DataFrame) -> dict:
    """For missense variants: share of diagnostic and of background variants reaching each PP3 level."""
    mis = cases[cases["class"] == "missense"]
    out = {}
    for score, cuts in PP3.items():
        if score not in mis:
            continue
        row = {}
        for label, sub in (("diagnostic", mis[mis["diagnostic"] == 1]), ("background", mis[mis["diagnostic"] == 0])):
            v = sub[score]
            row[label] = {"n": int(len(v)), **{lvl: round(float((v >= c).mean()), 4)
                                               for lvl, c in zip(["supporting", "moderate", "strong"], cuts)}}
        out[score] = row
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--study", default="snv")
    p.add_argument("--gnomad-max-af", type=float)
    p.add_argument("--grpmax", action="store_true")
    a = p.parse_args()
    data = DATA if a.study == "snv" else DATA / a.study
    cases, singles = load(data, evo2=a.study != "replication")
    before = (cases["case"].nunique(), len(cases))
    cases = gnomad_filter(cases, a.gnomad_max_af, a.grpmax)
    if a.study == "indel":
        # Rank each diagnostic indel among the background indels in the same genome.
        cases = cases[~cases["is_snv"]].reset_index(drop=True)
    labels, groups = cases["diagnostic"].to_numpy(), cases["sample"].to_numpy()
    lookups = [s for s in singles if s != "Evo 2 7B" and s not in FM_TABLES]
    combos = {"Lookups combined": lookups}
    for extra in [s for s in singles if s in FM_TABLES] + (["Evo 2 7B"] if "Evo 2 7B" in singles else []):
        combos[f"Lookups + {extra.replace(' 7B', '')}"] = lookups + [extra]
    for name, cols in combos.items():
        cases[name] = cross_validated(features(cases, cols), labels, groups)
    scores = singles + list(combos)
    fixed = {}
    if a.study == "replication":
        train = gnomad_filter(load(DATA, evo2=False)[0], a.gnomad_max_af, a.grpmax)
        fixed = {f"{name} (locked)": cols for name, cols in combos.items()}
        for name, cols in fixed.items():
            cases[name] = locked(train, cases, cols)
        scores += list(fixed)

    diag = cases[cases["diagnostic"] == 1].set_index("case")
    rng = np.random.default_rng(SEED)
    comp = {s: competitors(cases, s) for s in scores}
    strata = {"all": diag.index} | {g: diag.index[diag["group"] == g] for g in sorted(diag["group"].unique())}
    result = {"cases": len(diag), "median_variants_per_case": float(cases.groupby("case").size().median()),
              "coverage_all_variants": {s: int(cases[s].notna().sum()) for s in singles},
              "coverage_diagnostic": {s: int(diag[s].notna().sum()) for s in singles},
              "pp3_missense": pp3(cases), "strata": {}}
    if a.gnomad_max_af is not None:
        result["gnomad_filter"] = {"field": "AF_grpmax" if a.grpmax else "AF", "max_af": a.gnomad_max_af,
                                   "cases_before": before[0], "variants_before": before[1],
                                   "variants_after": len(cases)}
    for name, ids in strata.items():
        result["strata"][name] = {"cases": len(ids),
                                  "scores": {s: summarise(comp[s].loc[ids], diag.loc[ids, "sample"], rng) for s in scores},
                                  "added_to_lookups": {c: paired(comp[c].loc[ids], comp["Lookups combined"], diag.loc[ids, "sample"], rng)
                                                       for c in combos if c != "Lookups combined"}}
        if fixed:
            base = "Lookups combined (locked)"
            result["strata"][name]["added_to_locked_lookups"] = {
                c: paired(comp[c].loc[ids], comp[base], diag.loc[ids, "sample"], rng) for c in fixed if c != base}
    suffix = "" if a.study == "snv" else f"_{a.study}"
    if a.gnomad_max_af is not None:
        suffix += f"_gnomad{a.gnomad_max_af:g}" + ("_grpmax" if a.grpmax else "")
    out = HERE / "results" / f"results{suffix}.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print("PP3 on missense:", json.dumps(result["pp3_missense"]))
    for name, s in result["strata"].items():
        print(f"\n{name} ({s['cases']} cases)")
        for score, r in s["scores"].items():
            print(f"  {score:22s} top1 {r['top_1']:.2f}  top5 {r['top_5']:.2f}  top10 {r['top_10']:.2f} "
                  f"{r['top_10_ci95']}  median competitors {r['median_competitors']:.0f}")
        for c, d in (s["added_to_lookups"] | s.get("added_to_locked_lookups", {})).items():
            print(f"  {c} minus lookups: " + "  ".join(f"top{k} {d[f'top_{k}']:+.3f} {d[f'top_{k}_ci95']}" for k in TOP))


if __name__ == "__main__":
    main()
