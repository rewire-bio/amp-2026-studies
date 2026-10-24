"""Consequence class for every indel in a study, from MANE Select transcripts, independent of any score table.

CADD's precomputed indel table only covers indels seen in gnomAD, so its consequence is missing more often for
diagnostic indels than for background ones. Annotating all indels here gives the combined models the same
information for both. Classes, most severe first across overlapping MANE Select transcripts:

  splice      overlaps the two intronic bases at an exon-intron boundary
  frameshift  overlaps coding sequence, length change not a multiple of 3
  inframe     overlaps coding sequence, length change a multiple of 3
  utr         overlaps an exon outside coding sequence
  intronic    inside the transcript
  other

  python rare-disease/annotate_indels.py indel     # data/indel/case_variants.tsv -> consequence_indels.tsv
"""
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
ORDER = ["splice", "frameshift", "inframe", "utr", "intronic", "other"]


def transcripts() -> dict[str, list[dict]]:
    tx = {}
    with gzip.open(DATA / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz", "rt") as h:
        for line in h:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS") or 'tag "MANE_Select"' not in f[8]:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            t = tx.setdefault(tid, {"chrom": f[0], "exons": [], "cds": []})
            t["exons" if f[2] == "exon" else "cds"].append((int(f[3]), int(f[4])))
    by_chrom = defaultdict(list)
    for t in tx.values():
        t["exons"].sort()
        t["start"], t["end"] = t["exons"][0][0], t["exons"][-1][1]
        by_chrom[t["chrom"]].append(t)
    return by_chrom


def overlaps(a: int, b: int, intervals) -> bool:
    return any(s <= b and a <= e for s, e in intervals)


def classify(t: dict, a: int, b: int, change: int) -> str:
    if b < t["start"] or a > t["end"]:
        return "other"
    exons = t["exons"]
    splice = [(e + 1, e + 2) for _, e in exons[:-1]] + [(s - 2, s - 1) for s, _ in exons[1:]]
    if overlaps(a, b, splice):
        return "splice"
    if overlaps(a, b, t["cds"]):
        return "frameshift" if change % 3 else "inframe"
    if overlaps(a, b, exons):
        return "utr"
    return "intronic"


def main(study: str) -> None:
    out_dir = DATA / study
    by_chrom = transcripts()
    rows = []
    with (out_dir / "case_variants.tsv").open() as h:
        for r in csv.DictReader(h, delimiter="\t"):
            ref, alt, pos = r["ref"], r["alt"], int(r["pos"])
            if len(ref) == len(alt) == 1:
                continue
            # Bases touched after the VCF anchor; an insertion sits between pos and pos + 1.
            a, b = (pos + 1, pos + len(ref) - 1) if len(ref) > 1 else (pos, pos + 1)
            classes = [classify(t, a, b, len(alt) - len(ref)) for t in by_chrom.get(r["chrom"], [])]
            rows.append([r["chrom"], pos, ref, alt, min(classes or ["other"], key=ORDER.index)])
    with (out_dir / "consequence_indels.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "indel_class"])
        w.writerows(rows)
    counts = defaultdict(int)
    for *_, c in rows:
        counts[c] += 1
    print(f"{study}: {len(rows)} indels annotated", dict(counts))


if __name__ == "__main__":
    main(sys.argv[1])
