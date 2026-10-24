"""Lift the BRCA1 SGE variants from GRCh37 to hg38 and check every reference base against hg38.

  uv run --with pyliftover python brca1-case-ranking/liftover_brca1.py --study ../evo2-cloud-operations \
      --chain hg19ToHg38.over.chain.gz --out brca1-case-ranking/brca1_hg38.tsv

The hg38 reference base is fetched from UCSC's DAS-free sequence API for chr17 positions in one region request.
Writes id (GRCh37 key), chrom, pos38, ref, alt for every variant whose base matches; reports any that do not.
"""
import argparse
import csv
import json
import urllib.request
from pathlib import Path

from pyliftover import LiftOver


def hg38_sequence(start: int, end: int) -> str:
    """chr17:start-end (1-based, inclusive) from the UCSC REST API."""
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr17;start={start - 1};end={end}"
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)["dna"].upper()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--study", type=Path, required=True)
    p.add_argument("--chain", required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    lo = LiftOver(a.chain)
    with (a.study / "data/derived/brca1_variants.tsv").open() as h:
        variants = list(csv.DictReader(h, delimiter="\t"))
    lifted = []
    for v in variants:
        hit = lo.convert_coordinate("chr17", int(v["pos"]) - 1)
        if not hit or hit[0][0] != "chr17" or hit[0][2] != "+":
            print("no liftover:", v["id"], hit)
            continue
        lifted.append((v, hit[0][1] + 1))
    lo_pos, hi_pos = min(p for _, p in lifted), max(p for _, p in lifted)
    seq = hg38_sequence(lo_pos, hi_pos)
    ok, bad = [], []
    for v, pos38 in lifted:
        base = seq[pos38 - lo_pos]
        (ok if base == v["ref"].upper() else bad).append((v, pos38, base))
    with a.out.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["id", "chrom", "pos38", "ref", "alt"])
        for v, pos38, _ in ok:
            w.writerow([v["id"], "chr17", pos38, v["ref"], v["alt"]])
    print(f"{len(variants)} variants, {len(lifted)} lifted, {len(ok)} reference bases match hg38, {len(bad)} do not")
    for v, pos38, base in bad[:10]:
        print("mismatch:", v["id"], pos38, "hg38 has", base)


if __name__ == "__main__":
    main()
