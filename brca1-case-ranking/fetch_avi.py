"""Fetch AlphaGenome Atlas AVI scores for the BRCA1 benchmark variants (hg38).

  uv run --with alphagenome python brca1-case-ranking/fetch_avi.py --list-scorers
  uv run --with alphagenome python brca1-case-ranking/fetch_avi.py --scorer <AVI scorer name> --out brca1-case-ranking/avi_brca1.tsv

Key: ALPHAGENOME_API_KEY or ~/.config/alphagenome/key. Atlas outputs are for non-commercial use only and must not
be used to train other models (AlphaGenome terms); the per-variant table is not committed.
"""
import argparse
import csv
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent


def client():
    from alphagenome.atlas import atlas

    key = os.environ.get("ALPHAGENOME_API_KEY") or (Path.home() / ".config/alphagenome/key").read_text().strip()
    return atlas.create(key)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--list-scorers", action="store_true")
    p.add_argument("--scorer")
    p.add_argument("--variants", type=Path, default=HERE / "brca1_hg38.tsv")
    p.add_argument("--out", type=Path, default=HERE / "avi_brca1.tsv")
    a = p.parse_args()
    c = client()
    if a.list_scorers:
        for m in c.scorer_metadata():
            print(m)
        return
    from alphagenome.data import genome

    with a.variants.open() as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    variants = [genome.Variant(chromosome=r["chrom"], position=int(r["pos38"]),
                               reference_bases=r["ref"], alternate_bases=r["alt"], name=r["id"]) for r in rows]
    result = c.query_variants(variants, requested_scorers=[a.scorer])
    with a.out.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["id", "scorer", "value"])
        for scorer, ad in result.items():
            obs = ad.obs
            values = ad.X.reshape(len(obs), -1)
            for i, name in enumerate(obs.index if "name" not in obs else obs["name"]):
                w.writerow([name, scorer, repr(float(values[i, 0]))])
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
