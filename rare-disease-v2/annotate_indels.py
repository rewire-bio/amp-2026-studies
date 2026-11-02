"""Consequence class for indels, from MANE Select transcripts, independent of any score table.

Implements the indel column of PROTOCOL.md section 4. CADD's precomputed indel table only covers indels seen in
gnomAD, so its consequence would be missing more often for diagnostic indels than for background ones. Annotating
all indels here gives every indel the same information. Classes, most severe first across overlapping MANE Select
transcripts:

  splice      overlaps the two intronic bases at an exon-intron boundary
  frameshift  overlaps coding sequence, length change not a multiple of 3
  inframe     overlaps coding sequence, length change a multiple of 3
  utr         overlaps an exon outside coding sequence
  intronic    inside the transcript
  other

SNV rows in the input are skipped. The input needs chrom, pos, ref and alt columns; the output has those plus
indel_class.

  cd rare-disease-v2
  uv run -q --python 3.12 --with-requirements requirements.txt python annotate_indels.py --study indel
      # data/indel/case_variants.tsv to data/indel/consequence_indels.tsv
  uv run -q --python 3.12 --with-requirements requirements.txt python annotate_indels.py \\
      --variants ../rare-disease/data/clinvar_indel_plp2.tsv --out data/pool/consequence_indels.tsv

Derived from rare-disease/annotate_indels.py at amp-2026-v4: same classes and rules; paths from common.py, and
--variants/--out so it also annotates the ClinVar indel pool.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import common

ORDER = ["splice", "frameshift", "inframe", "utr", "intronic", "other"]


def transcripts() -> dict[str, list[dict]]:
    tx = {}
    with gzip.open(common.MANE_GTF, "rt") as h:
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


def annotate(variants: Path, out: Path) -> dict[tuple, str]:
    """Write chrom, pos, ref, alt, indel_class for every indel in `variants`; return the classes by key."""
    by_chrom = transcripts()
    rows = []
    with variants.open() as h:
        for r in csv.DictReader(h, delimiter="\t"):
            ref, alt, pos = r["ref"], r["alt"], int(r["pos"])
            if len(ref) == len(alt) == 1:
                continue
            # Bases touched after the VCF anchor; an insertion sits between pos and pos + 1.
            a, b = (pos + 1, pos + len(ref) - 1) if len(ref) > 1 else (pos, pos + 1)
            classes = [classify(t, a, b, len(alt) - len(ref)) for t in by_chrom.get(r["chrom"], [])]
            rows.append([r["chrom"], pos, ref, alt, min(classes or ["other"], key=ORDER.index)])
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(["chrom", "pos", "ref", "alt", "indel_class"])
        w.writerows(rows)
    counts = defaultdict(int)
    for *_, c in rows:
        counts[c] += 1
    print(f"{variants}: {len(rows)} indels annotated", dict(counts))
    return {(c, p, r, a): k for c, p, r, a, k in rows}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--study", help="shortcut for data/<study>/case_variants.tsv to data/<study>/consequence_indels.tsv")
    ap.add_argument("--variants", type=Path, help="tab-separated input with chrom, pos, ref, alt")
    ap.add_argument("--out", type=Path, help="output table")
    a = ap.parse_args()
    if a.study:
        variants = a.variants or common.study_dir(a.study) / "case_variants.tsv"
        out = a.out or common.study_dir(a.study) / "consequence_indels.tsv"
    elif a.variants and a.out:
        variants, out = a.variants, a.out
    else:
        ap.error("give --study, or both --variants and --out")
    annotate(variants, out)


if __name__ == "__main__":
    main()
