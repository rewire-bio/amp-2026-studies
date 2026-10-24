#!/usr/bin/env python3
"""Rank a case's variants by each score and report coverage.

A missing score is written as "missing" and ranked last for that score. It is never treated as zero,
because a missing CADD or SpliceAI value means "not scored", not "benign".
"""
import argparse
import csv

SCORES = {"CADD_RAW": "cadd", "SPLICEAI_MAX": "spliceai", "EVO2_LOF": "evo2"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    with open(args.scores) as handle:
        header = [h.split("]")[-1].split(":")[0].replace("INFO/", "") for h in handle.readline().lstrip("#").split("\t")]
        rows = [dict(zip([h.strip() for h in header], line.rstrip("\n").split("\t"))) for line in handle]
    for row in rows:
        row["id"] = f"{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}"
        for key in SCORES:
            row[key] = None if row[key] in (".", "") else float(row[key])
    for key, name in SCORES.items():
        scored = sorted((r for r in rows if r[key] is not None), key=lambda r: -r[key])
        for rank, row in enumerate(scored, start=1):
            row[f"{name}_rank"] = rank
        for row in rows:
            row.setdefault(f"{name}_rank", "missing")
        print(f"{name}: {len(scored)} of {len(rows)} variants scored")
    rows.sort(key=lambda r: r["cadd_rank"] if r["cadd_rank"] != "missing" else len(rows) + 1)
    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["variant", "cadd_raw", "cadd_rank", "spliceai_max", "spliceai_rank", "evo2_lof", "evo2_rank"])
        for r in rows:
            writer.writerow([r["id"],
                             "missing" if r["CADD_RAW"] is None else r["CADD_RAW"], r["cadd_rank"],
                             "missing" if r["SPLICEAI_MAX"] is None else r["SPLICEAI_MAX"], r["spliceai_rank"],
                             "missing" if r["EVO2_LOF"] is None else f"{r['EVO2_LOF']:.6f}", r["evo2_rank"]])


if __name__ == "__main__":
    main()
