"""Build Evo 2 input windows for every rare-disease case variant (GRCh38, 8,192 bp, uppercase).

One reference window per position and one variant window per allele, centred as in Arc's BRCA1 notebook
(the variant sits at index 4,096). For an indel the variant window keeps the reference length: REF is replaced by
ALT and the right flank is refilled from the genome, so both windows are 8,192 bp and mean log-likelihoods compare
like for like. Windows are uppercased: the UCSC FASTA marks repeats in lowercase and
Evo 2's byte tokenizer treats a and A as different inputs.

  python rare-disease/evo2/make_windows.py [--variants case_variants.tsv] [--out data] [--fasta hg38.fa.gz]
                                 [--skip-scored rare-disease/evo2/run]   # leave out variants with a score already
"""
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE.parent / "data/case_variants.tsv"
OUT = HERE / "data"
HG38 = HERE.parents[2] / "evo2-size-by-use-case/data/raw/hg38.fa.gz"
WINDOW = 8192


def records(path: Path, wanted: set[str]):
    """Yield (name, sequence) for the wanted chromosomes, one at a time."""
    name, parts = None, []
    with gzip.open(path, "rt") as h:
        for line in h:
            if line.startswith(">"):
                if name in wanted:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            elif name in wanted:
                parts.append(line.strip())
    if name in wanted:
        yield name, "".join(parts)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--variants", type=Path, default=CASES)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--fasta", type=Path, default=HG38)
    p.add_argument("--skip-scored", type=Path, action="append", default=[])
    a = p.parse_args()
    scored = set()
    for d in a.skip_scored:
        for f in d.glob("scores*.tsv"):
            with f.open() as h:
                scored |= {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in csv.DictReader(h, delimiter="\t")}
    by_chrom = defaultdict(list)
    skipped = 0
    with a.variants.open() as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if (r["chrom"], int(r["pos"]), r["ref"], r["alt"]) in scored:
                skipped += 1
                continue
            by_chrom[r["chrom"]].append((int(r["pos"]), r["ref"], r["alt"]))
    print(f"{skipped} variants already scored, left out")
    a.out.mkdir(exist_ok=True)
    n_ref = n_var = 0
    with gzip.open(a.out / "refs.tsv.gz", "wt") as refs, gzip.open(a.out / "vars.tsv.gz", "wt") as vars_:
        refs.write("ref_id\tseq\n")
        vars_.write("chrom\tpos\tref\talt\tref_id\tseq\n")
        for chrom, seq in records(a.fasta, set(by_chrom)):
            for pos in sorted({p for p, _, _ in by_chrom[chrom]}):
                p = pos - 1
                start, end = max(0, p - WINDOW // 2), min(len(seq), p + WINDOW // 2)
                window = seq[start:end].upper()
                offset = min(WINDOW // 2, p)
                ref_id = f"{chrom}:{pos}"
                refs.write(f"{ref_id}\t{window}\n")
                n_ref += 1
                for vpos, ref, alt in by_chrom[chrom]:
                    if vpos != pos:
                        continue
                    if window[offset:offset + len(ref)] != ref.upper():
                        raise ValueError(f"{chrom}:{pos} genome has {window[offset:offset + len(ref)]}, table has {ref}")
                    right = seq[p + len(ref):p + len(ref) + len(window) - offset - len(alt)].upper()
                    var_seq = window[:offset] + alt + right
                    if len(var_seq) != len(window):
                        raise ValueError(f"{chrom}:{pos} {ref}>{alt} window length {len(var_seq)}, reference {len(window)}")
                    vars_.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{ref_id}\t{var_seq}\n")
                    n_var += 1
            print(f"{chrom}: done", flush=True)
    print(f"{n_ref} reference windows, {n_var} variant windows, {n_ref + n_var} in total")


if __name__ == "__main__":
    main()
