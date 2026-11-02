"""Paths, keys and shared rules for the second run of the rare-disease case tests (PROTOCOL.md).

Every script imports this module, so the paths to the first run's inputs, the variant key, the study layout and the
gene lookup are defined once.
"""
from __future__ import annotations

import bisect
import csv
import os
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
# RDV2_DATA and RDV2_RESULTS point the scripts at another copy, for tests.
DATA = Path(os.environ.get("RDV2_DATA", HERE / "data"))
RESULTS = Path(os.environ.get("RDV2_RESULTS", HERE / "results"))
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


# ---------------------------------------------------------------------------------------------------------------
# MANE Select transcripts

def mane_transcripts() -> set[str]:
    """Ensembl transcript IDs (version stripped) of the MANE v1.5 Select transcripts."""
    import gzip
    import re
    out = set()
    with gzip.open(MANE_GTF, "rt") as h:
        for line in h:
            if '\ttranscript\t' in line and 'tag "MANE_Select"' in line:
                out.add(re.search(r'transcript_id "([^".]+)', line).group(1))
    return out


# ---------------------------------------------------------------------------------------------------------------
# CADD v1.7 by remote tabix (derived from rare-disease/lookup_scores.py cadd_file, amp-2026-v4)

CADD_SNV = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz"
CADD_INDEL = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel_inclAnno.tsv.gz"
SEVERITY = {c: i for i, c in enumerate(["STOP_GAINED", "STOP_LOST", "CANONICAL_SPLICE", "NON_SYNONYMOUS", "SPLICE_SITE",
                                        "SYNONYMOUS", "5PRIME_UTR", "3PRIME_UTR", "NONCODING_CHANGE", "INTRONIC",
                                        "UPSTREAM", "DOWNSTREAM", "REGULATORY", "INTERGENIC"])}


def cadd_lookup(wanted: set[tuple], workdir: Path, threads: int = 3) -> dict[tuple, dict]:
    """CADD PHRED, consequence and SpliceAI max delta for each wanted key (chr-prefixed).

    Consequence rule (PROTOCOL.md section 4): the most severe consequence over the rows for MANE Select transcripts;
    if the variant has no MANE row, the most severe over all rows (`mane` is then False). SpliceAI is the largest of
    the four delta scores over all rows. Indels come from CADD's gnomAD v4.0 indel table; indels not in it are absent.
    """
    import subprocess
    import time
    from concurrent.futures import ThreadPoolExecutor
    mane = mane_transcripts()
    workdir.mkdir(parents=True, exist_ok=True)
    rows: dict[tuple, dict] = {}
    snvs = {k for k in wanted if is_snv(k[2], k[3])}
    for url, subset, tag in ((CADD_SNV, snvs, "snv"), (CADD_INDEL, wanted - snvs, "indel")):
        if not subset:
            continue
        header = subprocess.check_output(["tabix", "-H", url], text=True, cwd=workdir).splitlines()[-1].lstrip("#")
        col = {c: i for i, c in enumerate(header.split("\t"))}
        splice = [col[c] for c in ("SpliceAI-acc-gain", "SpliceAI-acc-loss", "SpliceAI-don-gain", "SpliceAI-don-loss")]

        def one_chrom(chrom: str, subset=subset, url=url, tag=tag, col=col, splice=splice) -> dict:
            positions = sorted({p for c, p, _, _ in subset if c == chrom})
            regions = workdir / f"cadd_regions_{tag}_{chrom}.tsv"
            regions.write_text("".join(f"{chrom[3:]}\t{p}\t{p}\n" for p in positions))
            for attempt in range(8):
                r = subprocess.run(["tabix", "-R", str(regions), url], capture_output=True, text=True, cwd=workdir)
                if r.returncode == 0:
                    break
                time.sleep(5 * (attempt + 1))
            else:
                raise RuntimeError(f"CADD tabix failed for {chrom}: {r.stderr[-300:]}")
            found: dict[tuple, dict] = {}
            for line in r.stdout.splitlines():
                f = line.split("\t")
                key = ("chr" + f[0], int(f[1]), f[2], f[3])
                if key not in subset:
                    continue
                cons, gene = f[col["Consequence"]], f[col["GeneName"]]
                is_mane = f[col["FeatureID"]].split(".")[0] in mane
                spl = [float(f[i]) for i in splice if f[i] != "NA"]
                d = found.setdefault(key, {"cadd_phred": f[col["PHRED"]], "any": None, "any_gene": "",
                                           "mane_cons": None, "mane_gene": "", "spliceai_max": None})
                if d["any"] is None or SEVERITY.get(cons, 99) < SEVERITY.get(d["any"], 99):
                    d["any"], d["any_gene"] = cons, gene
                if is_mane and (d["mane_cons"] is None or SEVERITY.get(cons, 99) < SEVERITY.get(d["mane_cons"], 99)):
                    d["mane_cons"], d["mane_gene"] = cons, gene
                if spl and (d["spliceai_max"] is None or max(spl) > d["spliceai_max"]):
                    d["spliceai_max"] = max(spl)
            return found

        with ThreadPoolExecutor(threads) as pool:
            for found in pool.map(one_chrom, sorted({c for c, *_ in subset})):
                for k, d in found.items():
                    rows[k] = {"cadd_phred": d["cadd_phred"],
                               "cadd_consequence": d["mane_cons"] or d["any"],
                               "cadd_consequence_any": d["any"],
                               "gene": d["mane_gene"] or d["any_gene"],
                               "mane": d["mane_cons"] is not None,
                               "spliceai_max": d["spliceai_max"]}
    return rows


# ---------------------------------------------------------------------------------------------------------------
# gnomAD v3.1.2 genome counts from MyVariant.info (derived from rare-disease/gnomad_af.py, amp-2026-v4)

def hgvs(chrom: str, pos: int, ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}{ref}>{alt}"
    if alt == ref[0] and len(ref) > 1:  # deletion
        s, e = pos + 1, pos + len(ref) - 1
        return f"{chrom}:g.{s}del" if s == e else f"{chrom}:g.{s}_{e}del"
    if ref == alt[0] and len(alt) > 1:  # insertion
        return f"{chrom}:g.{pos}_{pos + 1}ins{alt[1:]}"
    e = pos + len(ref) - 1
    return f"{chrom}:g.{pos}_{e}delins{alt}" if e > pos else f"{chrom}:g.{pos}delins{alt}"


COUNT_FIELDS = ["ac", "an"] + [f"{x}_{g}" for g in GROUPS for x in ("ac", "an")]


def myvariant_counts(keys: list[tuple], batch: int = 1000) -> dict[tuple, dict]:
    """Full gnomAD v3.1.2 genome AC and AN, overall and per grpmax group, for each key found.

    Keys not in the result were not found in gnomAD genomes. Each value also has `af` and `grpmax` (the highest of the
    five group AFs, 0 where AN is 0).
    """
    import json
    import time
    import urllib.parse
    import urllib.request
    url = "https://myvariant.info/v1/variant?assembly=hg38&fields=gnomad_genome.ac,gnomad_genome.an"
    by_id = {hgvs(*k): k for k in keys}
    ids = list(by_id)
    out: dict[tuple, dict] = {}
    for i in range(0, len(ids), batch):
        body = urllib.parse.urlencode({"ids": ",".join(ids[i:i + batch])}).encode()
        for attempt in range(6):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, data=body), timeout=180) as r:
                    hits = json.load(r)
                break
            except OSError:
                time.sleep(5 * (attempt + 1))
        else:
            raise SystemExit("MyVariant.info did not answer")
        for h in hits:
            g = h.get("gnomad_genome")
            if not g or h.get("notfound"):
                continue
            if isinstance(g, list):
                g = g[0]
            ac, an = g.get("ac") or {}, g.get("an") or {}
            if isinstance(ac, list):
                ac = ac[0]
            if isinstance(an, list):
                an = an[0]
            d = {"ac": int(ac.get("ac", 0)), "an": int(an.get("an", 0))}
            for grp in GROUPS:
                d[f"ac_{grp}"] = int(ac.get(f"ac_{grp}", 0))
                d[f"an_{grp}"] = int(an.get(f"an_{grp}", 0))
            d["af"] = d["ac"] / d["an"] if d["an"] else 0.0
            d["grpmax"] = max((d[f"ac_{grp}"] / d[f"an_{grp}"] if d[f"an_{grp}"] else 0.0) for grp in GROUPS)
            out[by_id[h["query"]]] = d
        print(f"  MyVariant: {min(i + batch, len(ids))} of {len(ids)} queried, {len(out)} found", flush=True)
    return out
