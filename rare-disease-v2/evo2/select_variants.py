"""Choose the case variants Evo 2 must score, and reuse the first run's scores for the rest (PROTOCOL.md section 12).

Reads every study's case_variants.tsv and the first run's Evo 2 scores (../rare-disease/evo2/run/scores*.tsv and
run3/scores*.tsv; the first file with a key wins). Writes:
  evo2/run1_v1/scores.tsv     first-run scores for the case variants they cover (same columns)
  evo2/check_keys.tsv         200 of those keys, drawn with seed 20261111, scored again to check reproducibility
  evo2/variants_to_score.tsv  the variants with no first-run score, plus the 200 check keys

  cd rare-disease-v2 && python3 evo2/select_variants.py
then: uv run ... python evo2/make_windows.py --fasta data/spliceai/hg38.fa --variants evo2/variants_to_score.tsv
"""
import csv
import glob
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
FIRST = V2.parent / "rare-disease" / "evo2"
COLS = ["chrom", "pos", "ref", "alt"]


def key(r: dict) -> tuple:
    return (r["chrom"], r["pos"], r["ref"], r["alt"])


def order(k: tuple) -> tuple:
    return (k[0], int(k[1]), k[2], k[3])


def main() -> None:
    need = set()
    for study in ("snv", "indel", "replication"):
        with open(V2 / "data" / study / "case_variants.tsv") as h:
            need |= {key(r) for r in csv.DictReader(h, delimiter="\t")}
    old: dict[tuple, dict] = {}
    for f in sorted(glob.glob(str(FIRST / "run" / "scores*.tsv")) + glob.glob(str(FIRST / "run3" / "scores*.tsv"))):
        with open(f) as h:
            for r in csv.DictReader(h, delimiter="\t"):
                if key(r) in need and key(r) not in old:
                    old[key(r)] = r
    new = need - set(old)
    check = random.Random(20261111).sample(sorted(old), 200)
    (HERE / "run1_v1").mkdir(exist_ok=True)
    header = list(next(iter(old.values())))
    with open(HERE / "run1_v1" / "scores.tsv", "w") as h:
        h.write("\t".join(header) + "\n")
        for k in sorted(old, key=order):
            h.write("\t".join(old[k][c] for c in header) + "\n")
    with open(HERE / "check_keys.tsv", "w") as h:
        h.write("\t".join(COLS) + "\n")
        for k in sorted(check, key=lambda k: (k[0], int(k[1]))):
            h.write("\t".join(k) + "\n")
    with open(HERE / "variants_to_score.tsv", "w") as h:
        h.write("\t".join(COLS) + "\n")
        for k in sorted(new | set(check), key=order):
            h.write("\t".join(k) + "\n")
    print(f"{len(need)} case variants: {len(old)} reused, {len(new)} to score, plus {len(check)} check keys")


if __name__ == "__main__":
    main()
