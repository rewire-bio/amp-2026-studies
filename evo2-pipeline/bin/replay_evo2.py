#!/usr/bin/env python3
"""Replay recorded Evo 2 scores for the windows, so the pipeline runs without a GPU.

The recorded table comes from a GPU run (id, ref_score, var_score, delta). Every requested variant must be
present: a variant without a recorded score stops the run instead of silently becoming a missing value.
"""
import argparse
import csv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows", required=True)
    parser.add_argument("--recorded", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    with open(args.recorded) as handle:
        recorded = {row["id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with open(args.windows) as handle, open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", "ref_ll", "var_ll", "evo2_lof"])
        for row in csv.DictReader(handle, delimiter="\t"):
            key = f"chr{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}"
            if key not in recorded:
                raise SystemExit(f"No recorded Evo 2 score for {key}")
            r = recorded[key]
            writer.writerow([row["chrom"], row["pos"], row["ref"], row["alt"], r["ref_score"], r["var_score"],
                             repr(-float(r["delta"]))])


if __name__ == "__main__":
    main()
