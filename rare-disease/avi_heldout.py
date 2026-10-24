"""Does AVI's lead depend on its training data? AVI's labels are gnomAD frequency classes on chromosomes 1, 4, 7, 8,
10, 13 and 15 (Atlas preprint, p. 21), and the 1000 Genomes backgrounds are in gnomAD, so some background variants may
be AVI training examples. Repeat the single-score ranking within each split: cases whose diagnostic variant is on a
training chromosome, ranked only among background variants on training chromosomes, and the same for the other
chromosomes. Background sizes differ between the splits, so compare scores within a split, not across splits.

  uv run --with pandas --with numpy --with scikit-learn python rare-disease/avi_heldout.py   # results/avi_heldout.json
"""
import json

import case_benchmark as cb

TRAIN = {f"chr{c}" for c in (1, 4, 7, 8, 10, 13, 15)}
SCORES = ["CADD v1.7", "SpliceAI", "REVEL", "AlphaMissense", "GPN-MSA", "GPN-Star-V", "GPN-Star-M", "AlphaGenome AVI"]


def main() -> None:
    result = {}
    for study in ("snv", "replication"):
        cases, _ = cb.load(cb.DATA if study == "snv" else cb.DATA / study, evo2=False)
        result[study] = {}
        for split, on_train in (("training_chromosomes", True), ("other_chromosomes", False)):
            sub = cases[cases["chrom"].isin(TRAIN) == on_train]
            sub = sub[sub["case"].isin(sub.loc[sub["diagnostic"] == 1, "case"])]
            comp = {s: cb.competitors(sub, s) for s in SCORES}
            result[study][split] = {"cases": int(sub["diagnostic"].sum()),
                                    "median_variants_per_case": float(sub.groupby("case").size().median()),
                                    "top_1": {s: round(float((c < 1).mean()), 4) for s, c in comp.items()},
                                    "top_5": {s: round(float((c < 5).mean()), 4) for s, c in comp.items()}}
            r = result[study][split]
            print(f"{study} {split} ({r['cases']} cases, median {r['median_variants_per_case']:.0f} variants): "
                  + ", ".join(f"{s} {r['top_5'][s]:.2f}" for s in SCORES))
    (cb.HERE / "results/avi_heldout.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
