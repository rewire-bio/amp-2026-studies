"""Rare-disease case tests, second run: where does each score rank the diagnostic variant in a real background?

Each case (build_cases.py) is one 1000 Genomes person's rare variants in the panel regions plus one ClinVar P/LP
diagnostic variant. For every score, count the background variants that score equal to or higher than the
diagnostic (ties count against the score). A missing score ranks lowest, as a failed lookup would. The rule-based
scores are not a total order (scores from different tools are never compared), so their competitors are counted by
comparing each background variant with its case's diagnostic (section 7, rule_outcomes). Top k is the share of
cases with fewer than k competitors.

Implements PROTOCOL.md sections 2 and 3 (frequency and background inheritance filters), 4 (model classes), 7
(baselines), 8 (models, contrasts and the reading rule), 9 (checks), 10 (EVEE), 11 (intervals) and 12 (shared
resamples, so adding a score never moves another row).

Analyses (--analysis):
  main     grpmax AF below 0.1% for every variant, then the background inheritance filter, which counts the
           person's SNVs and indels together; then only the study's own variant type is kept (every analysis)
  bg1pct   grpmax AF below 1%, then the inheritance filter
  nomoi    the main frequency filter, no inheritance filter
  all1kg   background frequency with all of 1000 Genomes left out (grpmax_no1kg) below 0.1%, then the inheritance
           filter; diagnostics keep their full gnomAD grpmax
  heldout  main filters; each case ranked only among background variants on chromosomes of its diagnostic's split
           (AVI training chromosomes or the others); snv and replication only
  evee     bg1pct filters; lists restricted to variants with an EVEE row, cases whose diagnostic has none dropped;
           descriptive; snv and replication only
  sweep    exploratory frequency sweep (section 13) on the 5%-cap build in data/af5/<study>: background and
           diagnostic grpmax below each of common.SWEEP_AF, then the inheritance filter; snv and replication only.
           Writes results/<study>_sweep.json and results/analysis_sweep.md (sweep() below)

Scores: the single scores; Lookups combined (logistic regression on CADD, SpliceAI, REVEL and AlphaMissense with
model class and score-by-class terms, as in the first run); Lookups + X for each added score; Rule-based and
Rule + AVI (ACMG points, section 7). snv and indel are fitted with 5-fold cross-validation grouped by person.
replication is fitted that way on the 500 and also locked: fitted on the 200 SNV cases under the main filters (the
1% filter for bg1pct and evee) and applied unchanged (section 8). For heldout and evee the models are fitted on the
full filtered lists; the lists are then restricted for ranking.

Every analysis but heldout and evee also has the strata "gnomad: found" and "gnomad: absent" (the diagnostic found
in full gnomAD or not, section 9) and a gnomad_split block: the primary and key secondary contrasts in each, on
shared draws, and the difference absent minus found, flagged "possibly_gnomad_presence" when a gain (primary or
key secondary) is larger in the absent cases and that interval excludes zero. Contrasts for the two strata are
given only there, so each has one interval. Found means in gnomAD with AC above 0.

Input format assumptions: the cases.tsv `gene` and `moi` columns hold one value per gene, comma-separated (a single
moi value applies to every gene); moi is one of mono, biallelic, xl_mono, xl_bi or unknown (moi.py); any other
value is kept by the filter and counted in the output.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python case_benchmark.py \\
      --study snv|indel|replication --analysis main|bg1pct|nomoi|all1kg|heldout|evee [--drop-score NAME]

Writes results/<study>_<analysis>.json (common.RESULTS; with --drop-score, <study>_<analysis>_without_<name>.json).

Derived from rare-disease/case_benchmark.py at amp-2026-v4.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from importlib import metadata

import numpy as np
import pandas as pd

import common
from common import KEY

STUDY_CHOICES = ["snv", "indel", "replication"]
ANALYSES = ["main", "bg1pct", "nomoi", "all1kg", "heldout", "evee", "sweep"]
SNV_ONLY_ANALYSES = {"heldout", "evee", "sweep"}
BOOTSTRAP = 1000
TOP = [1, 5, 10]
SMALL_STRATUM = 20      # groups with fewer cases are flagged as small (section 8)
MIN_CLASS_DIAGNOSTICS = 10
GAIN_OF_INTEREST = 0.05  # smallest gain of interest, 5 points (section 8)

# Score files: (scores_<source>.tsv, column, name used in results).
SOURCES = [("cadd", "cadd_phred", "CADD"), ("cadd", "spliceai_max", "SpliceAI"), ("revel", "revel", "REVEL"),
           ("am", "alphamissense", "AlphaMissense"), ("avi", "avi", "AVI"), ("gpnmsa", "gpn_msa", "GPN-MSA"),
           ("gpnstar_v", "gpnstar_v", "GPN-Star-V"), ("gpnstar_m", "gpnstar_m", "GPN-Star-M"),
           ("evee", "evee", "EVEE")]
EVO2 = "Evo 2"
LOOKUPS = ["CADD", "SpliceAI", "REVEL", "AlphaMissense"]
ADDED = ["AVI", "GPN-MSA", "GPN-Star-V", "GPN-Star-M", EVO2]
GPN_FAMILY = ["GPN-MSA", "GPN-Star-V", "GPN-Star-M"]  # secondary, Bonferroni-adjusted over three
RULE, RULE_AVI = "Rule-based", "Rule + AVI"
RULE_INPUTS = {RULE: ["CADD", "SpliceAI", "REVEL", "AlphaMissense"],
               RULE_AVI: ["CADD", "SpliceAI", "REVEL", "AlphaMissense", "AVI"]}
LOCKED = " (locked)"

# ---------------------------------------------------------------------------------------------------------------
# Consequence classes (section 4)

CADD_CLASS = {"STOP_GAINED": "nonsense", "NON_SYNONYMOUS": "missense", "CANONICAL_SPLICE": "canonical_splice",
              "SPLICE_SITE": "splice_region", "SYNONYMOUS": "synonymous", "5PRIME_UTR": "utr", "3PRIME_UTR": "utr",
              "INTRONIC": "intronic"}  # anything else, STOP_LOST included, is "other"
INDEL_CLASS = {"frameshift": "frameshift", "inframe": "inframe", "splice": "canonical_splice", "utr": "utr",
               "intronic": "intronic"}  # annotate_indels.py; anything else is "other"
CLASSES = ["nonsense", "frameshift", "inframe", "missense", "canonical_splice", "splice_region", "synonymous",
           "utr", "intronic", "other"]
STUDY_CLASSES = {"snv": ["nonsense", "missense", "canonical_splice", "splice_region", "synonymous", "utr", "intronic",
                         "other"],
                 "indel": ["frameshift", "inframe", "canonical_splice", "utr", "intronic", "other"]}
ALWAYS_MERGED = ["synonymous", "utr"]  # merged into "other" for the fitted models

# Lookup tables whose indel coverage depends on gnomAD: 1000 Genomes is part of gnomAD, so background indels are
# almost always present and ClinVar diagnostic indels mostly absent. They are unavailable for every indel
# (section 6); SpliceAI for indels is the score we ran (run_spliceai.sh).
GNOMAD_INDEL_TABLES = ["CADD", "SpliceAI"]


def evo2_files() -> list:
    return sorted(common.EVO2_DIR.glob("run*/scores*.tsv"))


def load(study: str, drop: set[str], folder=None) -> tuple[pd.DataFrame, list[str], dict]:
    """Case variants with their scores and consequence class; returns (cases, single score names, load notes).

    folder: the study's files, common.study_dir(study) by default (the sweep reads common.sweep_dir(study)).
    """
    d = folder or common.study_dir(study)
    cases = pd.read_csv(d / "cases.tsv", sep="\t", dtype={"gene": str, "moi": str, "group": str, "class": str,
                                                          "zygosity": str, "clinvar_id": str})
    for c in ("gene", "moi", "group", "class", "zygosity"):
        cases[c] = cases[c].fillna("")
    cases["is_snv"] = (cases["ref"].str.len() == 1) & (cases["alt"].str.len() == 1)
    notes: dict = {}

    cadd_path = d / "scores_cadd.tsv"
    if not cadd_path.exists():
        raise SystemExit(f"{cadd_path} is missing; the SNV classes come from it")
    singles = []
    for source, col, name in SOURCES:
        path = d / f"scores_{source}.tsv"
        if not path.exists() or name in drop:
            continue
        extra = ["cadd_consequence"] if source == "cadd" and "cadd_consequence" not in cases else []
        s = pd.read_csv(path, sep="\t").drop_duplicates(KEY)
        if source == "cadd" and "mane" in s:
            # Section 4: variants whose consequence came from non-MANE rows, logged.
            in_cases = s[KEY].merge(cases[KEY].drop_duplicates(), on=KEY)
            no_mane = s.merge(in_cases, on=KEY)["mane"].astype(str).isin(["0", "False", "0.0"])
            notes["case_variants_without_mane_cadd_row"] = int(no_mane.sum())
        s = s[KEY + [col] + extra].rename(columns={col: name})
        if source == "evee":
            s["has_evee"] = True
        cases = cases.merge(s, on=KEY, how="left")
        singles.append(name)
    if "cadd_consequence" not in cases:  # CADD itself dropped: still read the consequence
        s = pd.read_csv(cadd_path, sep="\t")[KEY + ["cadd_consequence"]].drop_duplicates(KEY)
        cases = cases.merge(s, on=KEY, how="left")
    if "has_evee" in cases:
        cases["has_evee"] = cases["has_evee"].fillna(False).astype(bool)

    files = evo2_files() if EVO2 not in drop else []
    if files:
        e = pd.concat([pd.read_csv(f, sep="\t") for f in files]).drop_duplicates(KEY)[KEY + ["evo2_delta"]]
        e[EVO2] = -e["evo2_delta"]  # higher means more damaging
        cases = cases.merge(e[KEY + [EVO2]], on=KEY, how="left")
        singles.append(EVO2)

    # Class at analysis time (section 4): CADD's consequence for SNVs, annotate_indels.py for indels.
    cls = cases["cadd_consequence"].map(CADD_CLASS).fillna("other")
    indel_file = d / "consequence_indels.tsv"
    if indel_file.exists():
        ind = cases[KEY].merge(pd.read_csv(indel_file, sep="\t").drop_duplicates(KEY), on=KEY, how="left")
        ind_cls = ind["indel_class"].map(INDEL_CLASS).fillna("other")
        cls = cls.where(cases["is_snv"].to_numpy(), ind_cls.to_numpy())
    elif (~cases["is_snv"]).any():
        cls = cls.where(cases["is_snv"], "other")
        if study == "indel":  # in the SNV studies indels only count toward the inheritance filter, then are dropped
            notes["indels_without_annotation"] = int((~cases["is_snv"]).sum())
    cases["var_class"] = cls.to_numpy()
    diag = cases["diagnostic"] == 1
    if (cases.loc[diag, "class"] != "").any():
        notes["diagnostic_class_differs_from_cases_tsv"] = int(
            (cases.loc[diag, "class"] != cases.loc[diag, "var_class"]).sum())

    # Indels: CADD is never used, and SpliceAI is the run score (missing where it was not run).
    indel = ~cases["is_snv"]
    if "CADD" in cases:
        cases.loc[indel, "CADD"] = np.nan
    cases["spliceai_is_run"] = False
    if "SpliceAI" in cases:
        cases.loc[indel, "SpliceAI"] = np.nan
        run = d / "scores_spliceai_run.tsv"
        if run.exists():
            r = cases[KEY].merge(pd.read_csv(run, sep="\t").drop_duplicates(KEY), on=KEY, how="left")["spliceai_run"]
            cases.loc[indel, "SpliceAI"] = r[indel].to_numpy()
            cases.loc[indel, "spliceai_is_run"] = True
    return cases, singles, notes


# ---------------------------------------------------------------------------------------------------------------
# Filters (sections 2 and 3)

FILTERS = {
    "main": {"max_af": common.MAIN_AF, "background_column": "grpmax", "inheritance": True},
    "bg1pct": {"max_af": common.LOOSE_AF, "background_column": "grpmax", "inheritance": True},
    "nomoi": {"max_af": common.MAIN_AF, "background_column": "grpmax", "inheritance": False},
    "all1kg": {"max_af": common.MAIN_AF, "background_column": "grpmax_no1kg", "inheritance": True},
    "heldout": {"max_af": common.MAIN_AF, "background_column": "grpmax", "inheritance": True},
    "evee": {"max_af": common.LOOSE_AF, "background_column": "grpmax", "inheritance": True},
}
FILTERS_TEXT = {"main": "grpmax below 0.1%, inheritance filter", "bg1pct": "grpmax below 1%, inheritance filter"}
RESTRICTED = {"biallelic"}          # needs two hits for everyone
RESTRICTED_XX = {"xl_bi"}           # needs two hits for XX people
KNOWN_MOI = {"mono", "biallelic", "xl_mono", "xl_bi", "unknown", ""}
TWO_COPIES = {"hom", "hemi"}


def split_list(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def inheritance_keep(cases: pd.DataFrame, sex: dict[str, str]) -> tuple[np.ndarray, dict]:
    """Which rows the background inheritance filter keeps (section 3), and counts for the output.

    A background variant is dropped only if every gene it falls in needs two hits for this person and it has
    neither: it is not homozygous (or hemizygous) and the person has fewer than two variants in that gene after
    the frequency filter, the diagnostic included. Diagnostics and variants with no gene are always kept.
    """
    rows = []
    for i, gene, moi in zip(range(len(cases)), cases["gene"], cases["moi"]):
        genes, mois = split_list(gene), split_list(moi)
        if len(mois) == 1 and len(genes) > 1:
            mois = mois * len(genes)
        if genes and len(mois) != len(genes):
            raise SystemExit(f"row {i}: gene '{gene}' and moi '{moi}' do not pair up")
        rows += [(i, g, m.lower()) for g, m in zip(genes, mois)]
    hits = pd.DataFrame(rows, columns=["row", "gene", "moi"])
    hits["case"] = cases["case"].to_numpy()[hits["row"]]
    hits["sex"] = cases["sample"].map(sex).to_numpy()[hits["row"]]
    if hits["sex"].isna().any():
        missing = sorted(set(cases["sample"]) - set(sex))
        raise SystemExit(f"no sex in people.tsv for {missing[:5]}")
    zyg = cases["zygosity"].str.lower().to_numpy()[hits["row"]]
    n_in_gene = hits.groupby(["case", "gene"])["row"].transform("nunique").to_numpy()
    needs_two = hits["moi"].isin(RESTRICTED) | (hits["moi"].isin(RESTRICTED_XX) & (hits["sex"] == "female"))
    allowed = ~needs_two.to_numpy() | np.isin(zyg, list(TWO_COPIES)) | (n_in_gene >= 2)
    keep = np.ones(len(cases), dtype=bool)
    has_gene = np.zeros(len(cases), dtype=bool)
    has_gene[hits["row"].to_numpy()] = True
    any_allowed = np.zeros(len(cases), dtype=bool)
    any_allowed[hits.loc[allowed, "row"].to_numpy()] = True
    keep = ~has_gene | any_allowed | (cases["diagnostic"].to_numpy() == 1)
    unknown = sorted(set(hits["moi"]) - KNOWN_MOI)
    return keep, {"unrecognised_moi_values_kept": unknown, "rows_without_gene": int((~has_gene).sum())}


def apply_filters(cases: pd.DataFrame, analysis: str, study: str, sex: dict[str, str],
                  spec: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """Frequency filter, then inheritance filter, then drop the other variant type; cases keep their diagnostic.

    The inheritance filter counts the person's SNVs and indels together (section 3); then the snv and replication
    studies keep SNVs only and the indel study indels only. spec, in the layout of FILTERS, replaces the analysis's
    filters (the sweep's thresholds).
    """
    f = spec or FILTERS[analysis]
    diag = cases["diagnostic"].to_numpy() == 1
    bg_af = cases[f["background_column"]].to_numpy(float)
    dx_af = cases["grpmax"].to_numpy(float)
    af = np.where(diag, dx_af, bg_af)
    if np.isnan(af).any():
        raise SystemExit(f"{int(np.isnan(af).sum())} rows have no frequency in cases.tsv")
    info = {"max_af": f["max_af"], "background_frequency": f["background_column"],
            "diagnostic_frequency": "grpmax", "inheritance_filter": f["inheritance"],
            "cases_before": int(cases["case"].nunique()), "variants_before": len(cases)}
    keep = af < f["max_af"]
    if "kgp_max_af" in f:
        # The sweep's 1000 Genomes cap on background variants (section 13); diagnostics are not carried, so pass.
        if "kgp_af" not in cases:
            raise SystemExit("cases.tsv has no kgp_af column; the sweep needs the 5%-cap build (section 13)")
        kgp = cases["kgp_af"].to_numpy(float)
        keep &= diag | (kgp < f["kgp_max_af"])
        info["kgp_max_af"] = f["kgp_max_af"]
    out = cases[keep]
    solved = out.loc[out["diagnostic"] == 1, "case"]
    info["cases_lost_diagnostic"] = int(info["cases_before"] - solved.nunique())
    if info["cases_lost_diagnostic"]:
        print(f"WARNING: {info['cases_lost_diagnostic']} diagnostics fail the frequency filter; those cases are "
              "dropped", file=sys.stderr)
    out = out[out["case"].isin(solved)].reset_index(drop=True)
    info["variants_after_frequency"] = len(out)
    if f["inheritance"]:
        keep, notes = inheritance_keep(out, sex)
        out = out[keep].reset_index(drop=True)
        info |= notes
        info["variants_after_inheritance"] = len(out)
    # Each diagnostic is ranked among the person's background variants of its own type only, after filtering; the
    # other type counts toward the two-variant total above (section 3). A cases.tsv without the other type passes.
    own_type = ~out["is_snv"] if study == "indel" else out["is_snv"]
    other = "SNVs" if study == "indel" else "indels"
    info[f"background_{other.lower()}_in_inheritance_count"] = int((~own_type).sum()) if f["inheritance"] else 0
    out = out[own_type].reset_index(drop=True)
    info["variants_after_dropping_other_type"] = len(out)
    info["description"] = (
        f"background {f['background_column']} below {f['max_af']:g}"
        + (f" and 1000 Genomes AF below {f['kgp_max_af']:g}" if "kgp_max_af" in f else "")
        + f", diagnostics full gnomAD grpmax below {f['max_af']:g}"
        + ("; background inheritance filter" if f["inheritance"] else "; no inheritance filter")
        + f"; background {other} dropped after filtering")
    return out, info


# ---------------------------------------------------------------------------------------------------------------
# Model classes and features (section 4; features as in the first run)

def model_class_map(train: pd.DataFrame, kind: str) -> tuple[dict[str, str], dict]:
    """Merge synonymous, UTR and other, then every class with fewer than 10 training diagnostics, into other.

    kind is snv or indel: only that study type's classes are listed in the merges and the model classes; a class
    of the other type, which the ranked lists do not hold, maps to other.
    """
    dx = train.loc[train["diagnostic"] == 1, "var_class"]
    merged_first = dx.where(~dx.isin(ALWAYS_MERGED), "other")
    counts = merged_first.value_counts()
    own = STUDY_CLASSES[kind]
    small = [c for c in own if c not in ALWAYS_MERGED + ["other"] and counts.get(c, 0) < MIN_CLASS_DIAGNOSTICS]
    merged = [c for c in ALWAYS_MERGED if c in own] + small
    mapping = {c: ("other" if c in merged or c not in own else c) for c in CLASSES}
    info = {"training_diagnostics_per_class": {c: int(n) for c, n in dx.value_counts().sort_index().items()},
            "merged_into_other": merged, "model_classes": sorted({mapping[c] for c in own})}
    return mapping, info


def features(cases: pd.DataFrame, names: list[str]) -> pd.DataFrame:
    """Scores (missing set to 0), a missing flag on SNVs only, model class, and score x class interactions.

    "other" is the reference class. Missing flags are SNV-only: for indels, whether a table has a score reflects
    the study design, not anything a pipeline could use fairly.
    """
    out = pd.DataFrame(index=cases.index)
    cls = cases["model_class"]
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
    for c in sorted(set(CLASSES) - {"other"}):
        if not (cls == c).any():
            continue
        is_c = (cls == c).astype(float)
        out[f"is {c}"] = is_c
        for n in names:
            out[f"{n} x {c}"] = out[n] * is_c
    return out


def classifier():
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced"))


def cross_validated(cases: pd.DataFrame, names: list[str]) -> np.ndarray:
    """5-fold cross-validation grouped by person: no case is scored by a model that saw that person."""
    from sklearn.model_selection import GroupKFold

    x = features(cases, names).to_numpy()
    labels, groups = cases["diagnostic"].to_numpy(), cases["sample"].to_numpy()
    out = np.full(len(x), np.nan)
    for train, test in GroupKFold(n_splits=5).split(x, labels, groups):
        model = classifier().fit(x[train], labels[train])
        out[test] = model.predict_proba(x[test])[:, 1]
    return out


def locked(train: pd.DataFrame, test: pd.DataFrame, names: list[str]) -> np.ndarray:
    """Fit on every training case, then score the test cases with that model unchanged.

    Feature columns follow the training cases: a missing flag or class the training cases lack is dropped.
    """
    train = train.assign(**{n: np.nan for n in names if n not in train})  # a score the training set lacks
    x = features(train, names)
    model = classifier().fit(x.to_numpy(), train["diagnostic"].to_numpy())
    return model.predict_proba(features(test, names).reindex(columns=x.columns, fill_value=0.0).to_numpy())[:, 1]


# ---------------------------------------------------------------------------------------------------------------
# Rule-based baseline (section 7)

REVEL_PATHOGENIC = [(0.932, 4), (0.879, 3), (0.773, 2), (0.644, 1)]   # at or above
REVEL_BENIGN = [(0.016, -4), (0.052, -3), (0.183, -2), (0.290, -1)]   # at or below; 0.003 merged into -4
AM_PATHOGENIC = [(0.990, 4), (0.972, 3), (0.906, 2), (0.792, 1)]
AM_BENIGN = [(0.070, -3), (0.099, -2), (0.169, -1)]
PVS1_CLASSES = {"nonsense", "frameshift", "canonical_splice"}


def level(values: np.ndarray, pathogenic: list, benign: list) -> np.ndarray:
    """Points from calibrated thresholds; NaN where the score is missing."""
    out = np.where(np.isnan(values), np.nan, 0.0)
    for t, p in reversed(pathogenic):   # lowest threshold first, so the highest reached wins
        out = np.where(values >= t, p, out)
    for t, p in reversed(benign):       # highest threshold first, so the lowest reached wins
        out = np.where(values <= t, p, out)
    return out


def spliceai_level(values: np.ndarray) -> np.ndarray:
    """Walker et al. 2023 at supporting weight: +1 at 0.2 or more, -1 at 0.1 or less."""
    return np.where(np.isnan(values), np.nan, np.where(values >= 0.2, 1.0, np.where(values <= 0.1, -1.0, 0.0)))


def column(cases: pd.DataFrame, name: str) -> np.ndarray:
    return cases[name].to_numpy(float) if name in cases else np.full(len(cases), np.nan)


RULE_SCORES = (RULE, RULE_AVI)
LOSES, TIED, BEATS = 0, 1, 2  # one background variant against its case's diagnostic


def rule_keys(cases: pd.DataFrame) -> dict[str, np.ndarray]:
    """ACMG points per variant (section 7), the tool that set the level, that tool's score, CADD and AVI.

    points: +8 for nonsense, frameshift and canonical splice (tool none); missense the higher of the REVEL level
    (AlphaMissense where REVEL is missing) and the SpliceAI level, the tool being SpliceAI only when its level is
    strictly higher or the missense tools have no score; inframe the higher of +2 and the SpliceAI level; everything
    else the SpliceAI level (tool SpliceAI). Points are NaN where the class needs a score and has none. The tool is
    set for SNVs only; indels are compared on our SpliceAI score.
    """
    cls = cases["var_class"].to_numpy()
    snv = cases["is_snv"].to_numpy(bool)
    revel, am, spl = column(cases, "REVEL"), column(cases, "AlphaMissense"), column(cases, "SpliceAI")
    use_am = np.isnan(revel)
    mis_score = np.where(use_am, am, revel)
    mis_level = np.where(use_am, level(am, AM_PATHOGENIC, AM_BENIGN), level(revel, REVEL_PATHOGENIC, REVEL_BENIGN))
    mis_tool = np.where(np.isnan(mis_level), "none", np.where(use_am, "am", "revel"))
    spl_level = spliceai_level(spl)
    with np.errstate(invalid="ignore"):
        spl_wins = ~np.isnan(spl_level) & (np.isnan(mis_level) | (spl_level > mis_level))
    points = spl_level.copy()
    tool = np.where(np.isnan(spl_level), "none", "spliceai").astype(object)
    score = spl.copy()
    missense = cls == "missense"
    points[missense] = np.where(spl_wins, spl_level, mis_level)[missense]
    tool[missense] = np.where(spl_wins, "spliceai", mis_tool)[missense]
    score[missense] = np.where(spl_wins, spl, mis_score)[missense]
    inframe = cls == "inframe"
    points[inframe] = np.fmax(2.0, spl_level)[inframe]
    pvs1 = np.isin(cls, list(PVS1_CLASSES))
    points[pvs1] = 8.0
    tool[pvs1 | ~snv] = "none"
    score[pvs1 | ~snv] = np.nan
    return {"points": points, "tool": tool, "tool_score": score, "cadd": column(cases, "CADD"),
            "spliceai": spl, "avi": column(cases, "AVI"), "is_snv": snv}


def compare_higher(bg: np.ndarray, dx: np.ndarray) -> np.ndarray:
    """Higher wins (BEATS), equal or missing on both TIED, scored beats missing."""
    b, d = ~np.isnan(bg), ~np.isnan(dx)
    with np.errstate(invalid="ignore"):
        return np.where(b & d, np.where(bg > dx, BEATS, np.where(bg == dx, TIED, LOSES)),
                        np.where(b, BEATS, np.where(d, LOSES, TIED)))


def rule_outcomes(cases: pd.DataFrame, with_avi: bool) -> np.ndarray:
    """Each background variant against its case's diagnostic (section 7): LOSES, TIED or BEATS.

    The rule is not a total order (scores from different tools are never compared), so it is applied pair by pair.
    A background variant competes (TIED or BEATS) if it has more points, or equal points and:
      SNVs, same tool: a higher tool score, or an equal score and CADD PHRED at least as high;
      SNVs, different tools or no tool: CADD PHRED at least as high;
      a missing CADD on either side counts against the score (BEATS, not a tie);
      indels: our SpliceAI higher (BEATS), equal or missing on both (TIED); scored beats missing.
    TIED means it competes only because of equality. A diagnostic with no points ranks below every scored variant,
    and every background variant counts against it (BEATS), as a missing score does everywhere. Rule + AVI: at equal
    points a higher AVI wins (scored beats missing); equal or both missing AVI falls to the rule above.
    Diagnostic rows are LOSES. Each case must hold exactly one diagnostic.
    """
    k = rule_keys(cases)
    diag = cases["diagnostic"].to_numpy() == 1
    case = cases["case"].to_numpy()
    di = pd.Series(np.flatnonzero(diag), index=case[diag]).reindex(case)
    if di.isna().any():
        raise SystemExit("rule_outcomes: a case has no diagnostic")
    di = di.to_numpy(int)

    def dx(name: str) -> np.ndarray:
        return k[name][di]

    by_points = compare_higher(k["points"], dx("points"))  # TIED here means equal points (or none on both)
    c, e = k["cadd"], dx("cadd")
    with np.errstate(invalid="ignore"):
        by_cadd = np.where(np.isnan(c) | np.isnan(e), BEATS, np.where(c > e, BEATS, np.where(c == e, TIED, LOSES)))
        same_tool = (k["tool"] == dx("tool")) & (k["tool"] != "none")
        s, r = k["tool_score"], dx("tool_score")
        snv_eq = np.where(same_tool, np.where(s > r, BEATS, np.where(s == r, by_cadd, LOSES)), by_cadd)
    at_equal = np.where(dx("is_snv"), snv_eq, compare_higher(k["spliceai"], dx("spliceai")))
    if with_avi:
        by_avi = compare_higher(k["avi"], dx("avi"))
        at_equal = np.where(by_avi == TIED, at_equal, by_avi)
    out = np.where(by_points == TIED, at_equal, by_points)
    out = np.where(np.isnan(dx("points")), BEATS, out)
    return np.where(diag, LOSES, out)


def outcome_column(name: str) -> str:
    return f"{name} outcome"


# ---------------------------------------------------------------------------------------------------------------
# Ranks, resamples and intervals (section 11)

def competitors(cases: pd.DataFrame, score: str) -> pd.Series:
    """Background variants that count against the diagnostic, per case.

    A single score or model: those scoring at or above the diagnostic; all of them if the diagnostic has none. The
    rule scores: those whose stored pairwise outcome (rule_outcomes) is TIED or BEATS, so a subset of the rows (a
    heldout split, the EVEE lists) is counted from the same comparisons.
    """
    diag = cases["diagnostic"].to_numpy() == 1
    case = cases["case"].to_numpy()
    if score in RULE_SCORES:
        beats = cases[outcome_column(score)].to_numpy() >= TIED
        return pd.Series(beats & ~diag).groupby(case).sum().astype(int)
    s = cases[score].to_numpy(float)
    d = pd.Series(s[diag], index=case[diag])
    with np.errstate(invalid="ignore"):
        beats = ~diag & (s >= d.reindex(case).to_numpy())
    comp = pd.Series(beats).groupby(case).sum()
    background = pd.Series(~diag).groupby(case).sum()
    comp[d.isna().reindex(comp.index).to_numpy()] = background[d.isna().reindex(comp.index).to_numpy()]
    return comp.astype(int)


def tied(cases: pd.DataFrame, score: str) -> pd.Series:
    """Rule scores: whether at least one background variant is exactly tied with the diagnostic, per case.

    Exactly tied: it competes only because of equality (outcome TIED).
    """
    diag = cases["diagnostic"].to_numpy() == 1
    t = cases[outcome_column(score)].to_numpy() == TIED
    return pd.Series(t & ~diag).groupby(cases["case"].to_numpy()).any()


class Draws:
    """1,000 resamples of people, drawn once per stratum before any score is evaluated (section 11)."""

    def __init__(self, study: str, analysis: str, stratum: str, people) -> None:
        self.label = f"{study}/{analysis}/{stratum}"
        self.people = np.array(sorted(set(people)))
        seed = int.from_bytes(hashlib.sha256(self.label.encode()).digest()[:8], "big")
        rng = np.random.default_rng(seed)
        self.idx = rng.choice(len(self.people), size=(BOOTSTRAP, len(self.people)), replace=True)

    def means(self, values: np.ndarray, sample: np.ndarray) -> np.ndarray:
        """Mean of a per-case value in each resample: every drawn person brings all their cases."""
        pos = {p: i for i, p in enumerate(self.people)}
        where = np.array([pos[s] for s in sample], dtype=int)
        sums = np.bincount(where, weights=values, minlength=len(self.people))
        counts = np.bincount(where, minlength=len(self.people)).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            return sums[self.idx].sum(axis=1) / counts[self.idx].sum(axis=1)


def interval(boot: np.ndarray, level: float = 95.0) -> list[float]:
    tail = (100.0 - level) / 2
    return [round(float(x), 4) for x in np.nanpercentile(boot, [tail, 100.0 - tail])]


def summarise(comp: pd.Series, sample: pd.Series, draws: Draws) -> dict:
    c, s = comp.to_numpy(), sample.loc[comp.index].to_numpy()
    out = {"median_competitors": float(np.median(c))}
    for k in TOP:
        hit = (c < k).astype(float)
        out[f"top_{k}"] = round(float(hit.mean()), 4)
        out[f"top_{k}_ci95"] = interval(draws.means(hit, s))
    return out


BONFERRONI = round(100 - 5 / 3, 2)  # 98.33
BONFERRONI_EVO2 = round(100 - 5 / 4, 2)  # 98.75: Evo 2 was added later, so its row is adjusted over all four


def p_value(boot: np.ndarray) -> float:
    """Two-sided bootstrap p-value for a difference of zero: twice the smaller share of draws at or beyond zero.

    Floored at one draw's share, so it is never reported as 0. A difference that is zero in every draw gives 1.
    """
    b = boot[~np.isnan(boot)]
    if len(b) == 0:
        return float("nan")
    p = 2 * min(float(np.mean(b <= 0)), float(np.mean(b >= 0)))
    return round(min(1.0, max(p, 1 / len(b))), 4)


def paired(comp_a: pd.Series, comp_b: pd.Series, sample: pd.Series, draws: Draws, adjusted=False) -> dict:
    """Difference in top-k rate, a minus b, on the same cases, with the stratum's shared resamples.

    adjusted: False, True (the GPN family's Bonferroni level) or a level in percent (Evo 2's, over four scores).
    """
    a, b = comp_a.to_numpy(), comp_b.loc[comp_a.index].to_numpy()
    s = sample.loc[comp_a.index].to_numpy()
    out = {}
    for k in TOP:
        d = (a < k).astype(float) - (b < k).astype(float)
        boot = draws.means(d, s)
        out[f"top_{k}"] = round(float(d.mean()), 4)
        out[f"top_{k}_ci95"] = interval(boot)
        out[f"top_{k}_p"] = p_value(boot)
        if adjusted:
            level = BONFERRONI if adjusted is True else adjusted
            out[f"top_{k}_ci{level:g}"] = interval(boot, level)
    return out


def verdict(primary: dict, key_secondary: dict | None) -> str:
    """Section 8 reading rule for the primary contrast (top 1)."""
    est, (lo, hi) = primary["top_1"], primary["top_1_ci95"]
    same_way = key_secondary is not None and key_secondary["top_1"] > 0
    if lo > 0 and est >= GAIN_OF_INTEREST and same_way:
        return "AVI adds a gain that matters"
    if lo > -GAIN_OF_INTEREST and hi < GAIN_OF_INTEREST:
        return "No gain that matters"
    return "Inconclusive"


# ---------------------------------------------------------------------------------------------------------------
# Contrasts (section 8)

def contrast_specs(scores: list[str], study: str) -> list[dict]:
    """Every paired contrast available among the scores, with its role."""
    specs = []
    fits = [""] + ([LOCKED] if study == "replication" else [])
    for fit in fits:
        lookups, plus_avi = f"Lookups combined{fit}", f"Lookups + AVI{fit}"
        confirmatory = fit == LOCKED
        role = "primary" if confirmatory else "primary contrast, cross-validated (descriptive)"
        specs.append({"a": plus_avi, "b": lookups, "role": role, "primary": True, "confirmatory": confirmatory})
        specs.append({"a": lookups, "b": RULE, "role": "descriptive"})
        specs.append({"a": plus_avi, "b": RULE, "role": "descriptive"})
        for x in ADDED:
            if x == "AVI":
                continue
            gpn = x in GPN_FAMILY
            if x == "Evo 2":  # added after the other results; exploratory, adjusted over the four added scores
                specs.append({"a": f"Lookups + {x}{fit}", "b": lookups, "adjusted": BONFERRONI_EVO2,
                              "role": "exploratory, Bonferroni-adjusted over four scores"})
                continue
            specs.append({"a": f"Lookups + {x}{fit}", "b": lookups, "adjusted": gpn,
                          "role": ("secondary, Bonferroni-adjusted" if confirmatory else
                                   "secondary, cross-validated (descriptive)") if gpn else "descriptive"})
    specs.append({"a": RULE_AVI, "b": RULE, "role": "key secondary", "key_secondary": True})
    return [s for s in specs if s["a"] in scores and s["b"] in scores]


def evaluate_stratum(ids, diag: pd.DataFrame, comp: dict, specs: list[dict], draws: Draws,
                     with_contrasts: bool = True, ties: dict | None = None, with_verdict: bool = False) -> dict:
    """Top k for every score, and the paired contrasts, on one stratum's cases with its shared resamples.

    The verdict (section 8 reading rule) is given in the all-cases stratum only.
    """
    sample = diag["sample"]
    out = {"cases": len(ids), "people": int(diag.loc[ids, "sample"].nunique()), "small": len(ids) < SMALL_STRATUM,
           "scores": {s: summarise(c.loc[ids], sample, draws) for s, c in comp.items()}}
    if ties:
        out["rule_ties"] = {n: round(float(t.loc[ids].mean()), 4) for n, t in ties.items()}
    if not with_contrasts:
        return out
    contrasts = {}
    for spec in specs:
        row = {"role": spec["role"]} | paired(comp[spec["a"]].loc[ids], comp[spec["b"]].loc[ids], sample, draws,
                                              adjusted=spec.get("adjusted", False))
        contrasts[f"{spec['a']} minus {spec['b']}"] = row
    key = next((contrasts[f"{s['a']} minus {s['b']}"] for s in specs if s.get("key_secondary")), None)
    for spec in specs:
        if spec.get("primary"):
            row = contrasts[f"{spec['a']} minus {spec['b']}"]
            row["confirmatory"] = spec["confirmatory"]
            if with_verdict and spec["confirmatory"]:  # the reading rule is for the locked primary only
                row["verdict"] = verdict(row, key)
    out["contrasts"] = contrasts
    return out


def strata_of(diag: pd.DataFrame) -> dict[str, pd.Index]:
    """All cases, each diagnostic group, and each class within a group that mixes classes (the non-coding group)."""
    strata = {"all": diag.index}
    for g in sorted(diag["group"].unique()):
        ids = diag.index[diag["group"] == g]
        strata[g] = ids
        classes = diag.loc[ids, "var_class"]
        if classes.nunique() > 1:
            for c in sorted(classes.unique()):
                strata[f"{g}: {c}"] = ids[classes == c]
    return strata


# ---------------------------------------------------------------------------------------------------------------

def fit_models(cases: pd.DataFrame, train: pd.DataFrame | None, singles: list[str]) -> dict[str, list[str]]:
    """Add Lookups combined and Lookups + X columns (cross-validated, and locked when train is given)."""
    lookups = [s for s in LOOKUPS if s in singles]
    combos = {"Lookups combined": lookups}
    for x in ADDED:
        if x in singles:
            combos[f"Lookups + {x}"] = lookups + [x]
    for name, cols in combos.items():
        cases[name] = cross_validated(cases, cols)
    if train is not None:
        for name, cols in combos.items():
            cases[name + LOCKED] = locked(train, cases, cols)
    return combos


def package_versions() -> dict:
    out = {"python": sys.version.split()[0]}
    for p in ("numpy", "pandas", "scikit-learn", "scipy", "joblib", "threadpoolctl"):
        out[p] = metadata.version(p)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--study", choices=STUDY_CHOICES, required=True)
    p.add_argument("--analysis", choices=ANALYSES, default="main")
    p.add_argument("--drop-score", action="append", default=[], metavar="NAME",
                   help="leave a single score, and every row involving it, out (invariance test)")
    a = p.parse_args()
    study, analysis, drop = a.study, a.analysis, set(a.drop_score)
    if analysis in SNV_ONLY_ANALYSES and study == "indel":
        p.error(f"--analysis {analysis} covers the SNV studies only (PROTOCOL.md sections 9 and 10)")

    people = pd.read_csv(common.DATA / "people.tsv", sep="\t")
    sex = dict(zip(people["sample"], people["sex"].str.lower()))
    if bad := set(sex.values()) - {"male", "female"}:
        raise SystemExit(f"people.tsv sex values not male or female: {bad}")

    if analysis == "sweep":
        sweep(study, drop, sex)
        return

    cases, singles, notes = load(study, drop)
    cases, filt = apply_filters(cases, analysis, study, sex)
    train, train_filters = None, None
    if study == "replication":
        # Locked models are fitted under the main filters, except for the 1% check and EVEE (section 8).
        train_filters = "bg1pct" if analysis in ("bg1pct", "evee") else "main"
        train, _, _ = load("snv", drop)
        train, _ = apply_filters(train, train_filters, "snv", sex)
    # Model classes come from the training set: the 200 SNV cases (snv, replication) or the 120 indel cases.
    mapping, merges = model_class_map(train if study == "replication" else cases,
                                      "indel" if study == "indel" else "snv")
    merges["training_set"] = "the 120 indel cases" if study == "indel" else "the 200 SNV cases"
    cases["model_class"] = cases["var_class"].map(mapping)
    if train is not None:
        train["model_class"] = train["var_class"].map(mapping)

    if analysis != "evee":
        singles = [s for s in singles if s != "EVEE"]
    no_coverage = [s for s in singles if cases[s].notna().sum() == 0]
    singles = [s for s in singles if s not in no_coverage]

    # The ranked lists, fixed before any score is evaluated, and their resamples.
    ranked = cases
    if analysis == "evee":
        if "has_evee" not in cases:
            raise SystemExit("scores_evee.tsv is missing")
        ranked = cases[cases["has_evee"]]
        ranked = ranked[ranked["case"].isin(ranked.loc[ranked["diagnostic"] == 1, "case"])]
    diag = ranked[ranked["diagnostic"] == 1].set_index("case").sort_index()
    found, found_notes = None, {}
    if analysis == "heldout":
        strata = {"splits": diag.index}
    else:
        strata = strata_of(diag)
        if analysis != "evee":
            # Section 9: diagnostic found in full gnomAD or absent. A property of the cases, so these strata exist
            # whether or not a score is dropped; the contrasts in them exist only where the scores do.
            found, found_notes = gnomad_status(diag)
            for name, flag in GNOMAD_STRATA.items():
                if (found == flag).any():
                    strata[name] = diag.index[(found == flag).to_numpy()]
    draws = {name: Draws(study, analysis, name, diag.loc[ids, "sample"]) for name, ids in strata.items()}

    # Scores.
    combos = fit_models(cases, train, [s for s in singles if s != "EVEE"])
    if analysis == "evee":
        combos = {"Lookups combined": combos["Lookups combined"]}
    rules = []
    if analysis != "evee":
        # Rule + AVI only where AVI has coverage (not in the indel study, where it would equal Rule-based).
        rules = [n for n in RULE_SCORES if not set(RULE_INPUTS[n]) & drop and (n != RULE_AVI or "AVI" in singles)]
        for n in rules:
            cases[outcome_column(n)] = rule_outcomes(cases, with_avi=n == RULE_AVI)
    scores = singles + list(rules) + list(combos) + ([c + LOCKED for c in combos] if train is not None else [])
    ranked = cases.loc[ranked.index] if analysis == "evee" else cases

    result = {"study": study, "analysis": analysis,
              "filter": filt, "load_notes": notes,
              "cases": len(diag), "people": int(diag["sample"].nunique()),
              "median_variants_per_case": float(ranked.groupby("case").size().median()),
              "coverage": {s: {"variants": int(ranked[s].notna().sum()), "diagnostics": int(diag[s].notna().sum())}
                           for s in singles},
              "scores_without_coverage": no_coverage, "dropped_scores": sorted(drop),
              "model_classes": merges,
              "fits": "cross-validated (5-fold by person)" + (
                  f"; locked: fitted on the 200 SNV cases under the {train_filters} filters "
                  f"({FILTERS_TEXT[train_filters]})" if train is not None else ""),
              "intervals": f"{BOOTSTRAP} resamples of people per stratum, seed sha256('{study}/{analysis}/<stratum>')",
              "package_versions": package_versions()}
    if analysis == "evee":
        result["evee_note"] = ("Descriptive only (PROTOCOL.md section 10): agreement with ClinVar labels on in-sample "
                               "variants, not diagnostic performance. Models fitted on the full bg1pct lists, then "
                               "ranked within the EVEE-covered lists.")
        result["cases_without_evee_diagnostic"] = int(cases["case"].nunique() - len(diag))

    if analysis == "heldout":
        result["heldout"] = heldout(ranked, scores, study, draws["splits"])
    else:
        comp = {s: competitors(ranked, s) for s in scores}
        diag_ranked = ranked[ranked["diagnostic"] == 1].set_index("case").sort_index()
        ties = {n: tied(ranked, n) for n in rules}
        result["rule_ties"] = {n: round(float(t.mean()), 4) for n, t in ties.items()}
        specs = contrast_specs(scores, study) if analysis != "evee" else []
        result["strata"] = {name: evaluate_stratum(ids, diag_ranked, comp, specs, draws[name],
                                                   analysis != "evee" and name not in GNOMAD_STRATA,
                                                   ties, with_verdict=name == "all")
                            for name, ids in strata.items()}
        for name, s in result["strata"].items():
            s["seed_label"] = draws[name].label
        if found is not None:
            result["gnomad_split"] = found_notes | gnomad_split(comp, diag_ranked, found, study, analysis)
        primary = [k for k, v in result["strata"]["all"].get("contrasts", {}).items() if v.get("confirmatory")]
        primary = primary or [k for k, v in result["strata"]["all"].get("contrasts", {}).items() if "verdict" in v]
        if primary:  # keyed by the contrast's name, so the row names the scores it involves
            result["primary"] = {primary[0]: result["strata"]["all"]["contrasts"][primary[0]]}

    suffix = "".join(f"_without_{re.sub('[^A-Za-z0-9]+', '_', d)}" for d in sorted(drop))
    common.RESULTS.mkdir(parents=True, exist_ok=True)
    out = common.RESULTS / f"{study}_{analysis}{suffix}.json"
    out.write_text(json.dumps(clean(result), indent=2) + "\n")
    report(result)
    print(f"wrote {out}")


GNOMAD_STRATA = {"gnomad: found": True, "gnomad: absent": False}


def gnomad_status(diag: pd.DataFrame) -> tuple[pd.Series, dict]:
    """Whether each diagnostic is found in full gnomAD (section 9), indexed by case; None where unknown.

    Found means the site is in gnomAD with AC above 0. From data/gnomad_full.tsv, then data/pool/gnomad_full.tsv
    for diagnostics the first table lacks.
    """
    found = pd.Series([None] * len(diag), index=diag.index, dtype=object)
    keys = diag[KEY].reset_index()
    used = []
    for path in (common.DATA / "gnomad_full.tsv", common.DATA / "pool" / "gnomad_full.tsv"):
        if not path.exists():
            continue
        t = pd.read_csv(path, sep="\t", usecols=KEY + ["found", "ac"], dtype={"found": str}).drop_duplicates(KEY)
        site = t["found"].str.strip().str.lower().map({"1": True, "1.0": True, "true": True,
                                                       "0": False, "0.0": False, "false": False})
        # A site with no allele (AC 0) has no carrier in gnomAD, so it counts as absent, as in gnomad_overlap.md.
        t["flag"] = site.where(site.isna(), site.astype(bool) & (pd.to_numeric(t["ac"], errors="coerce") > 0))
        m = keys.merge(t[KEY + ["flag"]], on=KEY, how="left").set_index("case")["flag"].reindex(diag.index)
        fill = found.isna() & m.notna()
        found[fill] = m[fill]
        used.append(f"{path.parent.name}/{path.name}")
    notes = {"sources": used, "diagnostics_found": int((found == True).sum()),
             "diagnostics_absent": int((found == False).sum()),
             "diagnostics_without_status": int(found.isna().sum())}
    return found, notes


def gnomad_split(comp: dict, diag: pd.DataFrame, found: pd.Series, study: str, analysis: str) -> dict:
    """Section 9: primary and key secondary in cases whose diagnostic is found in gnomAD and absent from it.

    People are drawn once for all cases (seed label '<study>/<analysis>/gnomad split') and both strata are computed
    on the same draws, so the difference absent minus found is paired by person. A draw in which a stratum has no
    cases gives no value for it and is left out of the interval, as in the heldout check.
    """
    fit = LOCKED if study == "replication" else ""
    pairs = (("primary", f"Lookups + AVI{fit}", f"Lookups combined{fit}"), ("key secondary", RULE_AVI, RULE))
    draws = Draws(study, analysis, "gnomad split", diag["sample"])
    sample = diag["sample"]
    out, gains = {"seed_label": draws.label}, {}
    for name, flag in GNOMAD_STRATA.items():
        ids = diag.index[(found == flag).to_numpy()]
        row = {"cases": len(ids), "people": int(sample.loc[ids].nunique()), "contrasts": {}}
        for role, x, y in pairs:
            if x in comp and y in comp and len(ids):
                row["contrasts"][f"{x} minus {y}"] = {"role": role} | paired(comp[x].loc[ids], comp[y].loc[ids],
                                                                            sample, draws)
                gains[(name, role)] = ((comp[x].loc[ids] < 1).astype(float) - (comp[y].loc[ids] < 1).astype(float))
        out[name] = row
    out["difference"] = {}
    for role, x, y in pairs:
        if ("gnomad: found", role) not in gains or ("gnomad: absent", role) not in gains:
            continue
        ga, gf = gains[("gnomad: absent", role)], gains[("gnomad: found", role)]
        boot = (draws.means(ga.to_numpy(), sample.loc[ga.index].to_numpy())
                - draws.means(gf.to_numpy(), sample.loc[gf.index].to_numpy()))
        diff = {"role": role, "top_1_absent_minus_found": round(float(ga.mean() - gf.mean()), 4),
                "top_1_ci95": interval(boot), "top_1_p": p_value(boot)}
        diff["possibly_gnomad_presence"] = bool(diff["top_1_absent_minus_found"] > 0 and diff["top_1_ci95"][0] > 0)
        out["difference"][f"{x} minus {y}"] = diff
    return out


def heldout(cases: pd.DataFrame, scores: list[str], study: str, draws: Draws) -> dict:
    """Section 9: split cases by the diagnostic's chromosome (AVI training or not), rank within the split."""
    fit = LOCKED if study == "replication" else ""
    primary = (f"Lookups + AVI{fit}", f"Lookups combined{fit}")
    key = (RULE_AVI, RULE)
    out, gains = {}, {}
    for split, on_train in (("training_chromosomes", True), ("other_chromosomes", False)):
        sub = cases[cases["chrom"].isin(common.TRAIN_CHROMS) == on_train]
        sub = sub[sub["case"].isin(sub.loc[sub["diagnostic"] == 1, "case"])]
        diag = sub[sub["diagnostic"] == 1].set_index("case").sort_index()
        comp = {s: competitors(sub, s) for s in scores}
        row = {"cases": len(diag), "median_variants_per_case": float(sub.groupby("case").size().median()),
               "scores": {s: summarise(comp[s], diag["sample"], draws) for s in scores}, "contrasts": {}}
        for name, (x, y) in (("primary", primary), ("key secondary", key)):
            if x in comp and y in comp:
                row["contrasts"][f"{x} minus {y}"] = {"role": name} | paired(comp[x], comp[y], diag["sample"], draws)
                gains[(split, name)] = (((comp[x] < 1).astype(float) - (comp[y] < 1).astype(float)), diag["sample"])
        out[split] = row
    out["difference"] = {}
    for name, (x, y) in (("primary", primary), ("key secondary", key)):
        if ("training_chromosomes", name) not in gains:
            continue
        (gt, st), (go, so) = gains[("training_chromosomes", name)], gains[("other_chromosomes", name)]
        boot = draws.means(gt.to_numpy(), st.to_numpy()) - draws.means(go.to_numpy(), so.to_numpy())
        diff = {"contrast": f"{x} minus {y}", "top_1_training_minus_other": round(float(gt.mean() - go.mean()), 4),
                "top_1_ci95": interval(boot), "top_1_p": p_value(boot)}
        if name == "primary":
            diff["possibly_inflated_by_training_overlap"] = bool(diff["top_1_training_minus_other"] > 0
                                                                 and diff["top_1_ci95"][0] > 0)
        out["difference"][name] = diff
    out["note"] = ("Each split's resamples are the same draws of people (seed label '<study>/heldout/splits'), so "
                   "the difference between splits is paired by person.")
    return out


# ---------------------------------------------------------------------------------------------------------------
# Frequency sweep (section 13, exploratory)

SWEEP_LABELS = {common.MAIN_AF: "main", common.LOOSE_AF: "bg1pct"}  # these rows reuse those analyses' draws
AF_BANDS = [("below 0.1%", 0.0, common.MAIN_AF), ("0.1 to 1%", common.MAIN_AF, common.LOOSE_AF),
            ("1 to 5%", common.LOOSE_AF, common.SWEEP_BUILD_AF)]
WHY_SCORES = ["AVI", "CADD"]
SWEEP_TABLE = ["AVI", "CADD", "SpliceAI", "GPN-Star-V", EVO2]
SWEEP_OTHERS = ["REVEL", "AlphaMissense", "GPN-MSA", "GPN-Star-M"]
# The parts of a sweep row that must equal the main or bg1pct result at 0.1% and 1% (section 13).
REPRODUCED = ["cases", "people", "median_variants_per_case", "coverage", "rule_ties", "strata", "gnomad_split"]
NOT_REPRODUCED = {"verdict"}  # the sweep reads no verdict: it is exploratory


def sweep_label(af: float) -> str:
    """Seed label of a threshold: main and bg1pct at 0.1% and 1%, so those rows reproduce them, else sweep<af>."""
    return SWEEP_LABELS.get(af, f"sweep{af:g}")


def sweep_filter(af: float) -> dict:
    """Grpmax below af for background and diagnostic, 1000 Genomes AF below the larger of af and 1%, inheritance on.

    The 1000 Genomes cap is the one a build at that threshold would use (the v6 build used 1% for both), so at 1%
    and below the lists are the v6 lists.
    """
    return {"max_af": af, "background_column": "grpmax", "inheritance": True,
            "kgp_max_af": max(af, common.KGP_MAX_AF)}


def competitors_by_band(cases: pd.DataFrame, score: str) -> pd.DataFrame:
    """Competitors (as competitors(): at or above the diagnostic, all if it has no score), per case and grpmax band."""
    diag = cases["diagnostic"].to_numpy() == 1
    case = cases["case"].to_numpy()
    s = cases[score].to_numpy(float)
    d = pd.Series(s[diag], index=case[diag]).reindex(case).to_numpy()
    with np.errstate(invalid="ignore"):
        beats = ~diag & ((s >= d) | np.isnan(d))
    af = cases["grpmax"].to_numpy(float)
    return pd.DataFrame({band: pd.Series(beats & (af >= lo) & (af < hi)).groupby(case).sum()
                         for band, lo, hi in AF_BANDS})


def why(cases: pd.DataFrame) -> dict:
    """Spearman correlation of each WHY_SCORES score with background grpmax, and competitors by grpmax band.

    The correlation is over the distinct background variants of the ranked lists (each person's variant once, as a
    person's cases share a background); variants without the score are left out. Competitors by band are means per
    case, next to the mean number of background variants in each band.
    """
    from scipy.stats import spearmanr

    bg = cases[cases["diagnostic"] == 0].drop_duplicates(["sample"] + KEY)
    case = cases["case"].to_numpy()
    af = cases["grpmax"].to_numpy(float)
    background = cases["diagnostic"].to_numpy() == 0
    out = {"background_variants": len(bg), "spearman_with_grpmax": {}, "background_by_band": {}, "above_by_band": {}}
    for band, lo, hi in AF_BANDS:
        n = pd.Series(background & (af >= lo) & (af < hi)).groupby(case).sum()
        out["background_by_band"][band] = round(float(n.mean()), 3)
    for s in WHY_SCORES:
        if s not in cases:
            continue
        x = bg[[s, "grpmax"]].dropna()
        rho = spearmanr(x[s], x["grpmax"]).statistic if len(x) > 2 else np.nan
        out["spearman_with_grpmax"][s] = {"rho": round(float(rho), 4), "variants": len(x)}
        by = competitors_by_band(cases, s)
        out["above_by_band"][s] = {b: round(float(by[b].mean()), 3) for b in by}
    return out


def sweep_row(study: str, af: float, cases_all: pd.DataFrame, singles_all: list[str], train_all, sex: dict,
              drop: set[str]) -> dict:
    """One threshold: single scores, both baselines and AVI added to each, the gnomAD split, and the why columns.

    As main() for that filter, with only Lookups combined and Lookups + AVI fitted (adding a model changes no other
    row, section 11). Locked models (replication) are fitted on the 200 under the same threshold. Evo 2 is left out
    above 1%, where its scores do not cover the lists (section 13).
    """
    label = sweep_label(af)
    spec = sweep_filter(af)
    cases, filt = apply_filters(cases_all, label, study, sex, spec)
    train = None
    if train_all is not None:
        train, train_filt = apply_filters(train_all, label, "snv", sex, spec)
        filt["training_cases"] = train_filt["cases_before"] - train_filt["cases_lost_diagnostic"]
    mapping, merges = model_class_map(train if train is not None else cases, "snv")
    merges["training_set"] = "the 200 SNV cases"
    cases["model_class"] = cases["var_class"].map(mapping)
    if train is not None:
        train["model_class"] = train["var_class"].map(mapping)
    singles = [s for s in singles_all if s != "EVEE" and (s != EVO2 or af <= common.LOOSE_AF)]
    no_coverage = [s for s in singles if cases[s].notna().sum() == 0]
    singles = [s for s in singles if s not in no_coverage]

    diag = cases[cases["diagnostic"] == 1].set_index("case").sort_index()
    strata = {"all": diag.index}
    found, found_notes = gnomad_status(diag)
    for name, flag in GNOMAD_STRATA.items():
        if (found == flag).any():
            strata[name] = diag.index[(found == flag).to_numpy()]
    draws = {name: Draws(study, label, name, diag.loc[ids, "sample"]) for name, ids in strata.items()}

    lookups = [s for s in LOOKUPS if s in singles]
    combos = {"Lookups combined": lookups}
    if "AVI" in singles:
        combos["Lookups + AVI"] = lookups + ["AVI"]
    for name, cols in combos.items():
        cases[name] = cross_validated(cases, cols)
        if train is not None:
            cases[name + LOCKED] = locked(train, cases, cols)
    rules = [n for n in RULE_SCORES if not set(RULE_INPUTS[n]) & drop and (n != RULE_AVI or "AVI" in singles)]
    for n in rules:
        cases[outcome_column(n)] = rule_outcomes(cases, with_avi=n == RULE_AVI)
    scores = singles + rules + list(combos) + ([c + LOCKED for c in combos] if train is not None else [])

    comp = {s: competitors(cases, s) for s in scores}
    diag = cases[cases["diagnostic"] == 1].set_index("case").sort_index()
    ties = {n: tied(cases, n) for n in rules}
    specs = contrast_specs(scores, study)
    row = {"max_af": af, "seed_label_analysis": label, "filter": filt,
           "cases": len(diag), "people": int(diag["sample"].nunique()),
           "median_variants_per_case": float(cases.groupby("case").size().median()),
           "coverage": {s: {"variants": int(cases[s].notna().sum()), "diagnostics": int(diag[s].notna().sum())}
                        for s in singles},
           "scores_without_coverage": no_coverage, "model_classes": merges,
           "rule_ties": {n: round(float(t.mean()), 4) for n, t in ties.items()},
           "strata": {name: evaluate_stratum(ids, diag, comp, specs, draws[name], name not in GNOMAD_STRATA, ties)
                      for name, ids in strata.items()}}
    for name, s in row["strata"].items():
        s["seed_label"] = draws[name].label
    row["gnomad_split"] = found_notes | gnomad_split(comp, diag, found, study, label)
    row["why"] = why(cases)
    return row


def leaves(x, path=()) -> dict[tuple, object]:
    """Leaves of a JSON value keyed by their path; lists are leaves (intervals)."""
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            out |= leaves(v, path + (k,))
        return out
    return {path: x}


def reproduction(row: dict, reference: dict) -> dict:
    """Compare a sweep row with the main or bg1pct result: every REPRODUCED leaf the sweep row has must be equal.

    The reference has more rows (diagnostic groups, the GPN and Evo 2 models); those the sweep does not compute are
    not compared. A sweep leaf the reference lacks is a problem.
    """
    mine = leaves(clean({k: row[k] for k in REPRODUCED if k in row}))
    ref = leaves(clean({k: reference[k] for k in REPRODUCED if k in reference}))
    mine = {p: v for p, v in mine.items() if not set(p) & NOT_REPRODUCED}
    differ = [f"{' / '.join(map(str, p))}: {v} vs {ref[p]}" for p, v in mine.items()
              if p in ref and json.dumps(v) != json.dumps(ref[p])]
    missing = [" / ".join(map(str, p)) for p in mine if p not in ref]
    return {"compared": len(mine) - len(missing), "differ": len(differ), "missing_in_reference": len(missing),
            "identical": not differ and not missing, "examples": (differ + missing)[:20]}


def sweep(study: str, drop: set[str], sex: dict) -> None:
    """Section 13: every threshold of common.SWEEP_AF on the 5%-cap build; results/<study>_sweep.json.

    The 0.1% and 1% rows are compared with results/<study>_main.json and <study>_bg1pct.json (the v6 results), and
    the comparison is stored in the result. Then results/analysis_sweep.md is written from every sweep result.
    """
    folder = common.sweep_dir(study)
    cases_all, singles_all, notes = load(study, drop, folder)
    bg = cases_all["diagnostic"] == 0
    if (cases_all.loc[bg, "grpmax"] >= common.SWEEP_BUILD_AF).any():
        raise SystemExit(f"{folder}/cases.tsv has background variants at or above {common.SWEEP_BUILD_AF:g}")
    train_all = load("snv", drop, common.sweep_dir("snv"))[0] if study == "replication" else None
    result = {"study": study, "analysis": "sweep", "exploratory": (
                  "Added on 6 November 2026 after every other result had been seen (PROTOCOL.md section 13); not "
                  "confirmatory. Evo 2 is reported at 1% and below only."),
              "build": str(folder.relative_to(common.DATA.parent)) if folder.is_relative_to(common.DATA.parent)
              else str(folder), "build_max_af": common.SWEEP_BUILD_AF, "thresholds": common.SWEEP_AF,
              "load_notes": notes, "dropped_scores": sorted(drop),
              "fits": "cross-validated (5-fold by person)" + (
                  "; locked: fitted on the 200 SNV cases under the same threshold" if study == "replication" else ""),
              "intervals": f"{BOOTSTRAP} resamples of people per stratum, seed sha256('{study}/<label>/<stratum>'), "
                           "label main at 0.1%, bg1pct at 1%, sweep<threshold> otherwise",
              "package_versions": package_versions(), "rows": {}, "reproduction": {}}
    for af in common.SWEEP_AF:
        row = sweep_row(study, af, cases_all, singles_all, train_all, sex, drop)
        result["rows"][f"{af:g}"] = row
        print(sweep_line(row), flush=True)
        if af in SWEEP_LABELS:
            path = common.RESULTS / f"{study}_{SWEEP_LABELS[af]}.json"
            if path.exists():
                check = {"reference": path.name} | reproduction(row, json.loads(path.read_text()))
            else:
                check = {"reference": path.name, "identical": None, "note": "reference result not found"}
            result["reproduction"][f"{af:g}"] = check
            print(f"  reproduces {path.name}: {check.get('identical')} ({check.get('compared', 0)} values, "
                  f"{check.get('differ', 0)} differ, {check.get('missing_in_reference', 0)} missing)")
    common.RESULTS.mkdir(parents=True, exist_ok=True)
    suffix = "".join(f"_without_{re.sub('[^A-Za-z0-9]+', '_', d)}" for d in sorted(drop))
    out = common.RESULTS / f"{study}_sweep{suffix}.json"
    out.write_text(json.dumps(clean(result), indent=2) + "\n")
    print(f"wrote {out}")
    if not drop:
        write_sweep_report()


def sweep_line(row: dict) -> str:
    s = row["strata"]["all"]["scores"]
    return (f"sweep {row['max_af']:g}: {row['cases']} cases, median list {row['median_variants_per_case']:g}; "
            + "; ".join(f"{n} {s[n]['top_1']:.3f}" for n in ["AVI", "CADD", RULE] if n in s))


# Report: results/analysis_sweep.md, from results/<study>_sweep.json

def fmt_pct(x: float, signed: bool = False) -> str:
    x = round(x * 100, 1) + 0.0
    return f"{x:+.1f}" if signed and x != 0 else f"{x:.1f}"


def fmt_ci(r: dict, key: str = "top_1", signed: bool = False) -> str:
    lo, hi = r[f"{key}_ci95"]
    return f"{fmt_pct(r[key], signed)} ({fmt_pct(lo, signed)} to {fmt_pct(hi, signed)})"


def af_text(af: float) -> str:
    return f"{af * 100:g}%"


def sweep_tables(res: dict) -> list[str]:
    study = res["study"]
    fit = LOCKED if study == "replication" else ""
    combined, plus_avi = f"Lookups combined{fit}", f"Lookups + AVI{fit}"
    gain_combined, gain_rule = f"{plus_avi} minus {combined}", f"{RULE_AVI} minus {RULE}"
    built = next(iter(res["rows"].values()))["filter"]["cases_before"]
    fits = "locked, fitted on the 200 under the same threshold" if fit else "cross-validated"
    lines = [f"## {'Replication' if study == 'replication' else 'SNV'} study, {built} cases built", "",
             f"Top 1, % (95% interval). Combined: lookups combined ({fits}). Gains in percentage points with paired "
             "95% intervals. Evo 2: not scored above 1%.", "",
             "| Threshold | Cases | Median list | " + " | ".join(SWEEP_TABLE) + " | Rule-based | Combined | "
             "Rule + AVI minus Rule | Combined + AVI minus combined |", "|---" * (len(SWEEP_TABLE) + 7) + "|"]
    for key, row in res["rows"].items():
        st = row["strata"]["all"]
        s, con = st["scores"], st.get("contrasts", {})
        lost = row["filter"]["cases_lost_diagnostic"]
        cells = [af_text(row["max_af"]), f"{row['cases']}" + (f" ({lost} lost)" if lost else ""),
                 f"{row['median_variants_per_case']:g}"]
        cells += [fmt_ci(s[x]) if x in s else "not scored" for x in SWEEP_TABLE]
        cells += [fmt_ci(s[x]) if x in s else "" for x in (RULE, combined)]
        cells += [fmt_ci(con[g], signed=True) + f", p = {con[g]['top_1_p']:g}" if g in con else ""
                  for g in (gain_rule, gain_combined)]
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "Other single scores, top 1 %:", "", "| Threshold | " + " | ".join(SWEEP_OTHERS) + " |",
              "|---" * (len(SWEEP_OTHERS) + 1) + "|"]
    for row in res["rows"].values():
        s = row["strata"]["all"]["scores"]
        lines.append(f"| {af_text(row['max_af'])} | " + " | ".join(fmt_ci(s[x]) if x in s else "" for x in SWEEP_OTHERS)
                     + " |")
    lines += ["", "AVI alone by whether the diagnostic is in full gnomAD (AC above 0), top 1 %, and the gains in each "
              "(paired; absent minus found on shared draws):", "",
              "| Threshold | Found: cases | Found: AVI | Absent: cases | Absent: AVI | Rule + AVI gain, found | "
              "Rule + AVI gain, absent | Absent minus found |", "|---" * 8 + "|"]
    for row in res["rows"].values():
        st, gs = row["strata"], row["gnomad_split"]
        f, a = st.get("gnomad: found"), st.get("gnomad: absent")
        g = lambda name: gs.get(name, {}).get("contrasts", {}).get(gain_rule)  # noqa: E731
        d = gs.get("difference", {}).get(gain_rule)
        lines.append(f"| {af_text(row['max_af'])} | {f['cases'] if f else 0} | "
                     f"{fmt_ci(f['scores']['AVI']) if f and 'AVI' in f['scores'] else ''} | {a['cases'] if a else 0} | "
                     f"{fmt_ci(a['scores']['AVI']) if a and 'AVI' in a['scores'] else ''} | "
                     f"{fmt_ci(g('gnomad: found'), signed=True) if g('gnomad: found') else ''} | "
                     f"{fmt_ci(g('gnomad: absent'), signed=True) if g('gnomad: absent') else ''} | "
                     + (f"{fmt_pct(d['top_1_absent_minus_found'], True)} ({fmt_pct(d['top_1_ci95'][0], True)} to "
                        f"{fmt_pct(d['top_1_ci95'][1], True)})" if d else "") + " |")
    bands = [b for b, _, _ in AF_BANDS]
    lines += ["", "Why: Spearman correlation of each score with background grpmax (distinct background SNVs of the "
              "ranked lists), and the mean number of background variants per case ranked at or above the spiked-in "
              "variant, by background grpmax band, next to the mean number in each band:", "",
              "| Threshold | Background SNVs | AVI rho | CADD rho | " + " | ".join(f"In list, {b}" for b in bands)
              + " | " + " | ".join(f"Above by AVI, {b}" for b in bands) + " | "
              + " | ".join(f"Above by CADD, {b}" for b in bands) + " |", "|---" * (4 + 3 * len(bands)) + "|"]
    for row in res["rows"].values():
        w = row["why"]
        rho = lambda s: f"{w['spearman_with_grpmax'][s]['rho']:.3f}" if s in w["spearman_with_grpmax"] else ""  # noqa: E731
        above = lambda s: [f"{w['above_by_band'].get(s, {}).get(b, float('nan')):.2f}" for b in bands]  # noqa: E731
        lines.append(f"| {af_text(row['max_af'])} | {w['background_variants']} | {rho('AVI')} | {rho('CADD')} | "
                     + " | ".join(f"{w['background_by_band'][b]:.2f}" for b in bands) + " | "
                     + " | ".join(above("AVI")) + " | " + " | ".join(above("CADD")) + " |")
    lines += ["", "Reproduction of the v6 results by the 0.1% and 1% rows (every value both compute: cases, list "
              "length, coverage, ties, every score and contrast in all cases and in the gnomAD strata, and the gnomAD "
              "split):", ""]
    for key, check in res["reproduction"].items():
        if check.get("identical") is None:
            lines.append(f"- {af_text(float(key))}: {check['reference']} not found, not checked.")
        else:
            lines.append(f"- {af_text(float(key))} against `{check['reference']}`: "
                         + ("identical" if check["identical"] else "NOT identical")
                         + f", {check['compared']} values compared, {check['differ']} differ, "
                         f"{check['missing_in_reference']} missing from the reference.")
            lines += [f"  - {e}" for e in check["examples"]]
    return lines + [""]


def write_sweep_report() -> None:
    results = [json.loads(p.read_text()) for s in ("replication", "snv")
               if (p := common.RESULTS / f"{s}_sweep.json").exists()]
    if not results:
        return
    lines = ["# Frequency sweep (exploratory)", "",
             "PROTOCOL.md section 13. Added on 6 November 2026 after every other result had been seen, so it is "
             "exploratory, not confirmatory. Made by `case_benchmark.py --study <study> --analysis sweep` from the "
             f"5%-cap build in `data/{common.SWEEP_DATA}/`: `results/<study>_sweep.json`.", "",
             "At each threshold: background grpmax (family left out) below it, diagnostics full gnomAD grpmax below "
             "it (a case whose diagnostic fails is dropped), 1000 Genomes AF below the larger of the threshold and 1%, "
             "then the inheritance filter. Intervals are 1,000 resamples of people; the 0.1% and 1% rows use the "
             "main and bg1pct draws, so they reproduce those results.", ""]
    for res in results:
        lines += sweep_tables(res)
    text = "\n".join(lines)
    (common.RESULTS / "analysis_sweep.md").write_text(text)
    print(f"wrote {common.RESULTS / 'analysis_sweep.md'}")


def clean(x):
    """JSON-safe: NaN to null, numpy numbers to Python numbers."""
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return None if np.isnan(x) else float(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def report(result: dict) -> None:
    print(f"{result['study']} {result['analysis']}: {result['cases']} cases, median "
          f"{result['median_variants_per_case']:.0f} variants per case; {result['filter']['description']}")
    print("model classes merged into other:", result["model_classes"]["merged_into_other"])
    if "rule_ties" in result:
        print("diagnostics tied with a background variant:", result["rule_ties"])
    for k, v in result.get("gnomad_split", {}).get("difference", {}).items():
        flag = "  [possibly gnomAD presence]" if v.get("possibly_gnomad_presence") else ""
        print(f"  gnomAD absent minus found, {k}: {v['top_1_absent_minus_found']:+.3f} {v['top_1_ci95']}{flag}")
    if "heldout" in result:
        for split in ("training_chromosomes", "other_chromosomes"):
            r = result["heldout"][split]
            print(f"  {split} ({r['cases']} cases): " + "; ".join(
                f"{k} {v['top_1']:+.3f} {v['top_1_ci95']}" for k, v in r["contrasts"].items()))
        for k, v in result["heldout"]["difference"].items():
            print(f"  difference ({k}): {v['top_1_training_minus_other']:+.3f} {v['top_1_ci95']}")
        return
    s = result["strata"]["all"]
    for score, r in s["scores"].items():
        print(f"  {score:30s} top1 {r['top_1']:.3f} {r['top_1_ci95']}  top5 {r['top_5']:.3f}  top10 {r['top_10']:.3f}")
    for name, d in s.get("contrasts", {}).items():
        extra = f"  [{d['verdict']}]" if "verdict" in d else ""
        print(f"  {name}: top1 {d['top_1']:+.3f} {d['top_1_ci95']} ({d['role']}){extra}")


if __name__ == "__main__":
    main()
