"""Paths, keys and shared rules for the second run of the rare-disease case tests (PROTOCOL.md).

Every script imports this module, so the paths to the first run's inputs, the variant key, the study layout and the
gene lookup are defined once.
"""
from __future__ import annotations

import bisect
import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
V1 = HERE.parent / "rare-disease"           # the first run, read only
V1_DATA = V1 / "data"
BLOG = HERE.parents[1]                      # sibling repositories
HG38_GZ = BLOG / "evo2-size-by-use-case/data/raw/hg38.fa.gz"
REVEL_ZIP = BLOG / "evo2-cloud-operations/data/raw/current-predictors/revel-v1.3_all_chromosomes.zip"
MANE_GTF = V1_DATA / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
PANELAPP = V1_DATA / "panelapp_285.json"
REGION_GENES = V1_DATA / "region_genes.tsv"  # chrom, start, end, MANE gene symbol (exons plus 20 bp)

KEY = ["chrom", "pos", "ref", "alt"]
GROUPS = ["afr", "amr", "eas", "nfe", "sas"]  # gnomAD genome grpmax groups
MAIN_AF = 0.001    # main filter: grpmax AF below 0.1% (PROTOCOL.md section 2)
LOOSE_AF = 0.01    # looser background filter, and the build-time background filter
KGP_MAX_AF = 0.01  # 1000 Genomes AF below 1%, as in the first run
MAX_INDEL = 51     # 50 bp change plus the VCF anchor base

# Studies: which people they use and where their files live.
STUDIES = {"snv": "main", "indel": "main", "replication": "rep"}
SEEDS = {"snv": 20261102, "indel": 20261103, "replication": 20261109, "rep_people": 20261110}
TRAIN_CHROMS = {f"chr{c}" for c in (1, 4, 7, 8, 10, 13, 15)}  # AVI training chromosomes
PAR = {"chrX": [(10_001, 2_781_479), (155_701_383, 156_030_895)]}


def study_dir(study: str) -> Path:
    return DATA / study


def set_dir(people_set: str) -> Path:
    return DATA / people_set


def is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1


def in_par(chrom: str, pos: int) -> bool:
    return any(a <= pos <= b for a, b in PAR.get(chrom, []))


def haploid_region(chrom: str, pos: int, sex: str) -> bool:
    """True where an XY person has one copy: chrX outside the pseudoautosomal regions."""
    return chrom == "chrX" and sex == "male" and not in_par(chrom, pos)


class GeneIndex:
    """Variant position to the MANE genes whose padded exons contain it (region_genes.tsv)."""

    def __init__(self, path: Path = REGION_GENES):
        self.by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        with path.open() as h:
            for chrom, start, end, gene in csv.reader(h, delimiter="\t"):
                self.by_chrom[chrom].append((int(start), int(end), gene))
        for v in self.by_chrom.values():
            v.sort()
        self.starts = {c: [s for s, _, _ in v] for c, v in self.by_chrom.items()}
        self.longest = {c: max(e - s for s, e, _ in v) for c, v in self.by_chrom.items()}

    def genes(self, chrom: str, pos: int, end: int | None = None) -> set[str]:
        end = pos if end is None else end
        rows = self.by_chrom.get(chrom, [])
        if not rows:
            return set()
        lo = bisect.bisect_left(self.starts[chrom], pos - self.longest[chrom])
        hi = bisect.bisect_right(self.starts[chrom], end)
        return {g for s, e, g in rows[lo:hi] if s <= end and pos <= e}
