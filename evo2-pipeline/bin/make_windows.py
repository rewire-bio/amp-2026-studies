#!/usr/bin/env python3
"""Build reference and variant windows for SNVs, as in Arc's BRCA1 notebook.

Writes one row per variant: chrom, pos, ref, alt, window id, the reference window and the variant window.
Variants at the same position share a reference window. The reference base is checked against the genome.
"""
import argparse
import csv
import gzip


def read_fasta(path: str) -> str:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">"))


def snv_windows(chrom: str, pos: int, ref: str, alt: str, window: int) -> tuple[str, str]:
    p = pos - 1
    start = max(0, p - window // 2)
    ref_seq = chrom[start:min(len(chrom), p + window // 2)]
    offset = min(window // 2, p)
    if ref_seq[offset].upper() != ref.upper():
        raise SystemExit(f"Reference mismatch at {pos}: genome has {ref_seq[offset]}, VCF has {ref}")
    return ref_seq, ref_seq[:offset] + alt + ref_seq[offset + 1:]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snvs", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--window", type=int, default=8192)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    genome = read_fasta(args.reference)
    with open(args.snvs) as handle, open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", "ref_window", "ref_seq", "var_seq"])
        for chrom, pos, ref, alt in csv.reader(handle, delimiter="\t"):
            ref_seq, var_seq = snv_windows(genome, int(pos), ref, alt, args.window)
            writer.writerow([chrom, pos, ref, alt, f"{chrom}:{pos}", ref_seq, var_seq])


if __name__ == "__main__":
    main()
