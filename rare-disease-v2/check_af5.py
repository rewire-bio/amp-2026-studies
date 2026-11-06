"""The 5%-cap build against the v6 build: below the 1% caps, every case must hold the same rows (PROTOCOL.md section 13).

For each study (snv and replication), keeps the rows of data/af5/<study>/cases.tsv that a v6 build would have made:
the diagnostic, and background variants with grpmax below 1% and 1000 Genomes AF (kgp_af) below 1%. Drops kgp_af and
compares them with data/<study>/cases.tsv as text, case by case and in order. Then compares the scores of every v6
case variant: each scores_<source>.tsv row in data/<study>/ against the same key in data/af5/<study>/.

A differing case is listed with the rows found only in one build. For each such row the other build's values for the
same person and key are shown where it has them, so the cause can be read: a changed frequency, a 1000 Genomes AF
on the other side of 1%, or a variant missing from one build's person_variants.tsv.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python check_af5.py \\
      [--study snv replication]

Writes results/af5_identity.md (common.RESULTS). Exit status 0 when every case and every score is identical.
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

import common
from common import KEY

SOURCES = ["cadd", "revel", "am", "avi", "gpnmsa", "gpnstar_v", "gpnstar_m"]
SHOWN = ["case", "sample", *KEY, "diagnostic", "af", "grpmax", "grpmax_no1kg"]


def read(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def v6_rows(af5: pd.DataFrame) -> pd.DataFrame:
    """The af5 rows a v6 build would have made, without kgp_af."""
    diag = af5["diagnostic"] == "1"
    kgp = pd.to_numeric(af5["kgp_af"], errors="coerce")
    keep = diag | ((af5["grpmax"].astype(float) < common.LOOSE_AF) & (kgp < common.KGP_MAX_AF))
    return af5[keep].drop(columns="kgp_af").reset_index(drop=True)


def compare_cases(v6: pd.DataFrame, af5: pd.DataFrame) -> tuple[int, list[str]]:
    """Cases whose rows are identical, and a description of every other case."""
    if "kgp_af" not in af5:
        raise SystemExit("the af5 cases.tsv has no kgp_af column: build it with build_cases.py --cap 0.05")
    kept = v6_rows(af5)
    lines, same = [], 0
    if list(kept.columns) != list(v6.columns):
        return 0, [f"columns differ: {list(kept.columns)} vs {list(v6.columns)}"]
    v6_by, kept_by = dict(tuple(v6.groupby("case", sort=False))), dict(tuple(kept.groupby("case", sort=False)))
    af5_all = af5.set_index(["sample", *KEY])
    v6_all = v6.set_index(["sample", *KEY])
    for case in sorted(set(v6_by) | set(kept_by)):
        a, b = v6_by.get(case), kept_by.get(case)
        if a is not None and b is not None and a.reset_index(drop=True).equals(b.reset_index(drop=True)):
            same += 1
            continue
        lines.append(f"- {case}: " + ("missing from af5" if b is None else "missing from v6" if a is None else
                                      f"{len(a)} rows in v6, {len(b)} below the 1% caps in af5"))
        if a is None or b is None:
            continue
        rows_a, rows_b = {tuple(r) for r in a.to_numpy()}, {tuple(r) for r in b.to_numpy()}
        if rows_a == rows_b:
            lines.append("  - same rows in another order")
        for side, rows, other, other_name in (("only in v6", rows_a - rows_b, af5_all, "af5"),
                                              ("only in af5", rows_b - rows_a, v6_all, "v6")):
            for r in sorted(rows):
                d = dict(zip(a.columns, r))
                there = other.loc[[(d["sample"], d["chrom"], d["pos"], d["ref"], d["alt"])]] \
                    if (d["sample"], d["chrom"], d["pos"], d["ref"], d["alt"]) in other.index else None
                seen = "not there" if there is None else "; ".join(
                    f"{c} {there.iloc[0][c]}" for c in ("af", "grpmax", "grpmax_no1kg", "kgp_af") if c in there)
                lines.append(f"  - {side}: " + " ".join(d[c] for c in SHOWN) + f" (in {other_name}: {seen})")
    return same, lines


def compare_scores(study: str, keys: pd.DataFrame) -> list[str]:
    lines = []
    for source in SOURCES:
        old, new = common.study_dir(study) / f"scores_{source}.tsv", common.sweep_dir(study) / f"scores_{source}.tsv"
        if not old.exists() or not new.exists():
            lines.append(f"- {source}: {'v6' if not old.exists() else 'af5'} file missing, not compared")
            continue
        a, b = read(old).merge(keys, on=KEY), read(new).merge(keys, on=KEY)
        if set(a.columns) != set(b.columns):
            lines.append(f"- {source}: columns differ, {sorted(a.columns)} vs {sorted(b.columns)}")
            continue
        m = a.merge(b, on=KEY, how="outer", suffixes=("_v6", "_af5"), indicator=True)
        cols = [c for c in a.columns if c not in KEY]
        only = (m["_merge"] != "both").sum()
        both = m[m["_merge"] == "both"]
        differ = sum(int((both[f"{c}_v6"] != both[f"{c}_af5"]).sum()) for c in cols)
        lines.append(f"- {source}: {len(both)} v6 case variants scored in both, {only} scored in one build only, "
                     f"{differ} values differ")
    return lines


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--study", nargs="+", default=["snv", "replication"], choices=["snv", "replication"])
    a = p.parse_args()
    md = ["# The 5%-cap build against v6", "",
          "Made by `check_af5.py`. Rows of `data/af5/<study>/cases.tsv` with grpmax and 1000 Genomes AF below 1% "
          "(and every diagnostic), compared as text with `data/<study>/cases.tsv`, case by case and in order; then "
          "the scores of every v6 case variant in both builds.", ""]
    ok = True
    for study in a.study:
        v6, af5 = read(common.study_dir(study) / "cases.tsv"), read(common.sweep_dir(study) / "cases.tsv")
        same, lines = compare_cases(v6, af5)
        scores = compare_scores(study, v6[KEY].drop_duplicates())
        n = v6["case"].nunique()
        bad_scores = [x for x in scores if not x.endswith(", 0 scored in one build only, 0 values differ")]
        ok &= same == n and not lines and not bad_scores
        md += [f"## {study}", "",
               f"{same} of {n} cases identical below the 1% caps. The af5 build has {len(af5)} rows, "
               f"{len(af5) - len(v6_rows(af5))} of them background variants at or above a 1% cap.", ""]
        md += lines + ([""] if lines else []) + ["Scores of the v6 case variants:", ""] + scores + [""]
    md.append("Result: " + ("every case and every score identical." if ok else "DIFFERENCES, listed above."))
    text = "\n".join(md) + "\n"
    print(text)
    common.RESULTS.mkdir(parents=True, exist_ok=True)
    (common.RESULTS / "af5_identity.md").write_text(text)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
