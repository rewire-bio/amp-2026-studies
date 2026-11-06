"""Each person's rare variants in the panel regions, with zygosity: data/<set>/person_variants.tsv.

Implements the background rule of PROTOCOL.md section 5 (1000 Genomes AF below 1% in the regions) and the zygosity
that section 3 needs. Variant rules are those of the first run's background(): TYPE snp or indel, INFO AF below 1%,
one ALT allele, ref and alt at most 51 bases, ACGT only, and carried (the GT has a 1). A record repeated by
overlapping fetch regions is kept once. Zygosity is hom for two alt alleles, hemi for a haploid call (chrX outside
the PARs in XY people), het otherwise.

Inputs: the main set's genotypes from ../rare-disease/data/kgp/*.bcf and the replication set's from
data/rep/kgp/*.bcf (both made by ../rare-disease/fetch_kgp.sh). Writes data/main/person_variants.tsv and
data/rep/person_variants.tsv with columns sample, chrom, pos, ref, alt, zygosity, kgp_af.

  uv run -q --python 3.12 --with-requirements requirements.txt python backgrounds.py [main] [rep] [--cap AF]

--cap replaces the 1% 1000 Genomes cap; the frequency sweep (PROTOCOL.md section 13) runs it with --cap 0.05 and
RDV2_DATA=data/af5, so the v6 files are not touched (AF_SWEEP_RUN.md). The comparison with the first run is made at
the 1% cap only.

Derived from rare-disease/build_cases.py (background) at amp-2026-v4.
"""
from __future__ import annotations

import csv
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from common import DATA, KGP_MAX_AF, MAX_INDEL, V1_DATA, is_snv

SETS = {"main": V1_DATA / "kgp", "rep": DATA / "rep" / "kgp"}
TYPES = '(TYPE="snp" || TYPE="indel")'


def zygosity(gt: str) -> str:
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) == 1:
        return "hemi"
    return "hom" if all(a == "1" for a in alleles) else "het"


def background(kgp: Path, cap: float = KGP_MAX_AF) -> list[tuple]:
    rows, seen = [], set()
    bcfs = sorted(kgp.glob("chr*.bcf"))
    assert len(bcfs) == 23, f"{kgp}: {len(bcfs)} chromosome files, expected 23"
    for bcf in bcfs:
        samples = subprocess.run(["bcftools", "query", "-l", str(bcf)], capture_output=True, text=True,
                                 check=True).stdout.split()
        out = subprocess.run(["bcftools", "query", "-i", f"{TYPES} && AF<{cap}",
                              "-f", "%CHROM\t%POS\t%REF\t%ALT\t%AF[\t%GT]\n", str(bcf)],
                             capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            chrom, pos, ref, alt, af, *gts = line.split("\t")
            if "," in alt or max(len(ref), len(alt)) > MAX_INDEL or not set(ref + alt) <= set("ACGT"):
                continue
            for s, gt in zip(samples, gts):
                if "1" in gt and (s, chrom, pos, ref, alt) not in seen:  # overlapping regions repeat records
                    seen.add((s, chrom, pos, ref, alt))
                    rows.append((s, chrom, int(pos), ref, alt, zygosity(gt), af))
    return rows


def main() -> None:
    args = sys.argv[1:]
    cap = KGP_MAX_AF
    if "--cap" in args:
        i = args.index("--cap")
        cap = float(args[i + 1])
        args = args[:i] + args[i + 2:]
    for set_name, kgp in SETS.items():
        if args and set_name not in args:
            continue
        rows = background(kgp, cap)
        expected = set((DATA / set_name / "samples50.txt").read_text().split()) if set_name == "rep" else \
            set((V1_DATA / "samples50.txt").read_text().split())
        assert {r[0] for r in rows} == expected, f"{set_name}: samples differ from samples50.txt"
        out = DATA / set_name / "person_variants.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as h:
            w = csv.writer(h, delimiter="\t", lineterminator="\n")
            w.writerow(["sample", "chrom", "pos", "ref", "alt", "zygosity", "kgp_af"])
            w.writerows(rows)

        snv = defaultdict(int)
        indel = defaultdict(int)
        for r in rows:
            (snv if is_snv(r[3], r[4]) else indel)[r[0]] += 1
        zyg = Counter(r[5] for r in rows)
        hemi_where = Counter(r[1] for r in rows if r[5] == "hemi")
        print(f"{set_name}: {len(rows)} person variants, {len({r[1:5] for r in rows})} distinct; per person median "
              f"SNVs {statistics.median(snv.values())} (range {min(snv.values())} to {max(snv.values())}), "
              f"indels {statistics.median(indel.values())}; zygosity {dict(zyg)}; hemi calls by chromosome "
              f"{dict(hemi_where)}")

        if set_name == "main" and cap == KGP_MAX_AF:
            compare_first_run(rows)


def compare_first_run(rows: list[tuple]) -> None:
    """Main set SNVs, less ClinVar P/LP keys, against each person's background in the first run's 200 SNV cases."""
    with (V1_DATA / "clinvar_plp2.tsv").open() as h:
        plp = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in csv.DictReader(h, delimiter="\t")}
    ours = defaultdict(set)
    for s, chrom, pos, ref, alt, _, _ in rows:
        if is_snv(ref, alt) and (chrom, pos, ref, alt) not in plp:
            ours[s].add((chrom, pos, ref, alt))
    first = defaultdict(lambda: defaultdict(set))
    with (V1_DATA / "cases.tsv").open() as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["diagnostic"] == "0":
                first[r["sample"]][r["case"]].add((r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
    same = sum(all(keys == ours[s] for keys in cases.values()) for s, cases in first.items())
    sizes = [len(k) for cases in first.values() for k in cases.values()]
    print(f"main vs first run: {same} of {len(first)} people have exactly the first run's background SNVs in every "
          f"case (first-run median background per case {statistics.median(sizes)}, ours per person "
          f"{statistics.median(len(v) for v in ours.values())})")


if __name__ == "__main__":
    main()
