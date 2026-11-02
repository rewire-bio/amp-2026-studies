"""gnomAD v3.1.2 frequencies as a patient's pipeline would see them: full, family left out, all 1000 Genomes left out.

Implements PROTOCOL.md section 1 (frequencies, formula, checks) and the "all of 1000 Genomes left out" check of
section 9. Steps, each run by its own subcommand (`all` runs them in order):

  full     Full gnomAD v3.1.2 genome AC and AN, overall and per group, from MyVariant.info (common.myvariant_counts)
           for every distinct key in data/<set>/person_variants.tsv, the first run's diagnostics and ClinVar P/LP
           pools, and data/pool/*_pool.tsv if present. Incremental: keys already in data/gnomad_full.tsv are not
           queried again.
  fetch    gnomAD's own genotype calls for the person-variant positions, from the HGDP + 1KG subset release
           (gnomad.genomes.v3.1.2.hgdp_tgp.chrN.vcf.bgz on Google Cloud Storage). One process per chromosome, 14 at
           a time. Cached as slim VCFs (GT, GQ, DP, AD of the released 1000 Genomes samples and of the samples that
           make up the subset file's own AC and AN) in data/gnomad_subset/dense/chrN.partK.vcf.gz; a rerun fetches
           only positions not fetched before.
  counts   Splits multiallelic rows (bcftools norm -m-), applies gnomAD's adj filter and counts alleles per row:
           data/gnomad_subset/counts_chrN.tsv.
  af       data/<set>/person_variant_af.tsv: af_fam, grpmax_fam, af_no1kg, grpmax_no1kg, in_subset.
  checks   results/gnomad_checks.md.
  overlap  results/gnomad_overlap.md.

Why the fetch does not use bcftools view -R: the subset file is dense (4,151 samples per row, about 35 kB compressed
per row) and its tabix index has one entry per 16 kb window (about 23 MB of file), so each region query reads on
average half a window. The fetch reads the index itself, finds each position by bisection on BGZF block boundaries
with HTTP range requests, and reads only the rows at that position.

Adj filter (protocol section 1): GQ at least 20; DP at least 10; for a heterozygote, alt allele depth over DP at
least 0.2 (gnomAD's get_adj_expr). A genotype failing it counts for neither AC nor AN. Haploid calls (chrX outside
the PARs in XY samples) count 1 allele. Deviation: the protocol allows DP 5 for haploid calls, but the subset file's
own AC and AN hold haploid calls to DP 10 (with 5, 4 of 1,592 chrX rows outside the PARs reproduce them; with 10, all
do), so haploid calls are held to DP 10 here. The checks recount both ways.

Formula (family left out), for group g in afr, amr, eas, nfe, sas and overall:
  AF = max(AC_g - AC_fam,g, 0) / (AN_g - AN_fam,g), 0 if the denominator is not positive.
Family members count in their gnomad_population_inference group; members in fin or oth change only the overall AF.
Where the subset file has no row, AC_fam is 0 and AN_fam is 2 per released member in the group (1 for an XY member
on chrX outside the PARs). All 1000 Genomes left out: the same with all 2,500 released 1000 Genomes samples; where
there is no row, AN_g is reduced by the group's share of 1000 Genomes samples, AN_g * (1 - n_1kg,g / N_g), with N_g
the number of gnomAD v3.1.2 genomes in the group. A key not found in gnomAD has 0 for every frequency.

  uv run -q --python 3.12 --with-requirements requirements.txt python gnomad_af.py all

Derived from rare-disease/gnomad_af.py at amp-2026-v4.
"""
from __future__ import annotations

import bisect
import csv
import gzip
import json
import random
import statistics
import struct
import subprocess
import sys
import threading
import time
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

from common import DATA, GROUPS, KEY, RESULTS, V1_DATA, haploid_region, in_par, myvariant_counts

CACHE = DATA / "gnomad_subset"
DENSE = CACHE / "dense"
META = CACHE / "gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz"
PEDIGREE = CACHE / "1kGP.3202_samples.pedigree_info.txt"
FULL = DATA / "gnomad_full.tsv"
SETS = ["main", "rep"]
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]
GCS = "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"
SUBSET_URL = GCS + "/gnomad.genomes.v3.1.2.hgdp_tgp.{chrom}.vcf.bgz"
SITES_URL = GCS + "/gnomad.genomes.v3.1.2.sites.{chrom}.vcf.bgz"
# gnomAD v3.1.2 genomes per group (gnomAD v3.1 release notes; 76,156 in all)
N_GNOMAD = {"afr": 20744, "amr": 7647, "eas": 2604, "nfe": 34029, "sas": 2419}
N_GNOMAD_ALL = 76156
KGP_GROUPS = GROUPS + ["other"]
HAPLOID_DP = 10  # protocol says 5; the subset file's own counts use 10 (deviation, see the module docstring)
FULL_COLS = ["found", "ac", "an"] + [f"{x}_{g}" for g in GROUPS for x in ("ac", "an")] + ["af", "grpmax"]
INFO_KEEP = ["AC", "AN", "gnomad_AC", "gnomad_AN"] + [f"gnomad_{x}_{g}" for g in GROUPS for x in ("AC", "AN")]
FMT_KEEP = ["GT", "GQ", "DP", "AD"]


def log(*a) -> None:
    print(*a, flush=True)


# ---------------------------------------------------------------------------------------------------------------
# Samples and keys

def samples() -> dict[str, dict]:
    """gnomAD ID to project, release status, membership of the subset file's own AC/AN, group and (1000 Genomes) sex."""
    csv.field_size_limit(10**9)
    with PEDIGREE.open() as h:
        sex = {r[0]: {"1": "male", "2": "female"}[r[3]] for r in (line.split() for line in list(h)[1:])}
    out = {}
    with gzip.open(META, "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            proj = json.loads(r["hgdp_tgp_meta"])["project"]
            rel = r["relatedness_inference"]
            related = json.loads(rel).get("related") if rel.startswith("{") else None
            pop_json = r["gnomad_population_inference"]
            pop = json.loads(pop_json).get("pop") if pop_json.startswith("{") else ""
            kgp_id = r["s"][:-1] if r["s"][-1] in "AB" else r["s"]
            out[r["s"]] = {"project": proj, "release": r["gnomad_release"] == "true",
                           "subset": r["high_quality"] == "true" and related is False, "pop": pop,
                           "sex": sex.get(kgp_id, ""), "kgp_release": proj == "1000 Genomes"
                           and r["gnomad_release"] == "true"}
    assert sum(s["kgp_release"] for s in out.values()) == 2500
    return out


def people() -> dict[str, dict]:
    with (DATA / "people.tsv").open() as h:
        return {r["sample"]: {**r, "family": [x for x in r["family"].split(",") if x]}
                for r in csv.DictReader(h, delimiter="\t")}


def read_keys(path: Path, where=None) -> set[tuple]:
    with path.open() as h:
        return {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in csv.DictReader(h, delimiter="\t")
                if where is None or where(r)}


def person_variants(set_name: str) -> list[dict]:
    with (DATA / set_name / "person_variants.tsv").open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def present_sets() -> list[str]:
    return [s for s in SETS if (DATA / s / "person_variants.tsv").exists()]


def trim(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Minimal representation after a multiallelic split: common suffix, then common prefix, removed."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


# ---------------------------------------------------------------------------------------------------------------
# Step 1: full counts

def step_full() -> None:
    keys = set()
    for s in SETS:
        if (DATA / s / "person_variants.tsv").exists():
            keys |= read_keys(DATA / s / "person_variants.tsv")
    for f in ("cases.tsv", "replication/cases.tsv"):
        keys |= read_keys(V1_DATA / f, lambda r: r["diagnostic"] == "1")
    for f in ("clinvar_plp2.tsv", "clinvar_indel_plp2.tsv"):
        keys |= read_keys(V1_DATA / f)
    for f in sorted((DATA / "pool").glob("*_pool.tsv")):
        keys |= read_keys(f)
    have: dict[tuple, list] = {}
    if FULL.exists():
        with FULL.open() as h:
            for r in csv.DictReader(h, delimiter="\t"):
                have[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = [r[c] for c in FULL_COLS]
    new = sorted(keys - set(have))
    log(f"full: {len(keys)} keys, {len(new)} to query")
    found = myvariant_counts(new) if new else {}
    for k in new:
        d = found.get(k)
        have[k] = [1] + [d[c] for c in FULL_COLS[1:]] if d else [0] + [0] * (len(FULL_COLS) - 1)
    order = {c: i for i, c in enumerate(CHROMS)}
    with FULL.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(KEY + FULL_COLS)
        for k in sorted(have, key=lambda k: (order.get(k[0], 99), k[1], k[2], k[3])):
            w.writerow(list(k) + have[k])
    n_found = sum(int(v[0]) for v in have.values())
    log(f"full: {len(have)} keys written, {n_found} found in gnomAD genomes")


def load_full() -> dict[tuple, dict]:
    out = {}
    with FULL.open() as h:
        for r in csv.DictReader(h, delimiter="\t"):
            d = {c: (float(r[c]) if c in ("af", "grpmax") else int(r[c])) for c in FULL_COLS}
            out[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = d
    return out


# ---------------------------------------------------------------------------------------------------------------
# Step 2: fetch rows from the dense subset file by bisection on BGZF blocks

MAGIC = b"\x1f\x8b\x08\x04"


class Remote:
    """HTTP range reads with retries, one session per thread."""

    def __init__(self, url: str):
        self.url = url
        self.local = threading.local()
        self.bytes = 0
        self.requests = 0
        self.lock = threading.Lock()

    def get(self, start: int, end: int) -> bytes:
        if not hasattr(self.local, "s"):
            self.local.s = requests.Session()
        for attempt in range(8):
            try:
                r = self.local.s.get(self.url, headers={"Range": f"bytes={start}-{end - 1}"}, timeout=120)
                if r.status_code == 206 and 0 < len(r.content) <= end - start:
                    with self.lock:
                        self.bytes += len(r.content)
                        self.requests += 1
                    return r.content
                if r.status_code == 416:
                    return b""
            except requests.RequestException:
                pass
            time.sleep(min(60, 2 ** attempt))
            self.local.s = requests.Session()
        raise RuntimeError(f"range read failed: {self.url} {start}-{end}")


def block_at(buf: bytes, i: int) -> int | None:
    """Size of the BGZF block starting at buf[i], or None if it is not a complete block header."""
    if buf[i:i + 4] != MAGIC or buf[i + 10:i + 12] != b"\x06\x00" or buf[i + 12:i + 14] != b"BC":
        return None
    return struct.unpack("<H", buf[i + 16:i + 18])[0] + 1


def inflate(buf: bytes, i: int, size: int) -> bytes:
    return zlib.decompress(buf[i + 18:i + size - 8], -15)


def read_tbi(path: Path) -> dict[str, list[int]]:
    """Linear index (virtual offsets, one per 16 kb window) per contig of a tabix index."""
    d = gzip.open(path).read()
    assert d[:4] == b"TBI\x01"
    n_ref, *_, l_nm = struct.unpack("<8i", d[4:36])
    names = d[36:36 + l_nm].split(b"\0")[:n_ref]
    p = 36 + l_nm
    out = {}
    for name in names:
        n_bin, = struct.unpack("<i", d[p:p + 4])
        p += 4
        for _ in range(n_bin):
            _, n_chunk = struct.unpack("<Ii", d[p:p + 8])
            p += 8 + 16 * n_chunk
        n_intv, = struct.unpack("<i", d[p:p + 4])
        p += 4
        out[name.decode()] = list(struct.unpack(f"<{n_intv}Q", d[p:p + 8 * n_intv]))
        p += 8 * n_intv
    return out


class Reader:
    """Sequential reader of decompressed lines from a block offset of a remote BGZF file.

    skip is the number of bytes to drop from the first block (a tabix virtual offset), or None to drop everything up
    to the first newline (a block found by bisection, which may start inside a row).
    """
    CHUNK = 256 * 1024

    def __init__(self, remote: Remote, offset: int, skip: int | None):
        self.remote, self.offset = remote, offset  # offset: next compressed byte to fetch
        self.raw = b""
        self.text = b""
        self.pre = b""
        self.mode, self.skip = ("skip", skip) if skip is not None else ("newline", 0)
        self.it = None

    def position(self) -> int:
        """Compressed offset of the next block not yet decompressed."""
        return self.offset - len(self.raw)

    def _more(self) -> bool:
        while True:
            if len(self.raw) >= 18:
                size = block_at(self.raw, 0)
                if size is None:
                    raise RuntimeError("lost BGZF block boundary")
                if len(self.raw) >= size:
                    data = inflate(self.raw, 0, size)
                    self.raw = self.raw[size:]
                    if self.mode == "skip":
                        data, self.mode = data[self.skip:], "done"
                    elif self.mode == "newline":
                        self.pre += data
                        nl = self.pre.find(b"\n")
                        if nl < 0:
                            continue
                        data, self.pre, self.mode = self.pre[nl + 1:], b"", "done"
                    self.text += data
                    return True
            chunk = self.remote.get(self.offset, self.offset + self.CHUNK)
            if not chunk:
                return False
            self.offset += len(chunk)
            self.raw += chunk

    def lines(self):
        while True:
            nl = self.text.find(b"\n")
            while nl < 0:
                if not self._more():
                    return
                nl = self.text.find(b"\n")
            line, self.text = self.text[:nl], self.text[nl + 1:]
            yield line


def probe(remote: Remote, off: int) -> tuple[int, int] | None:
    """(compressed offset of the first block at or after off, POS of the first line that starts after a newline in it).

    Returns None if no such line exists before the end of the file.
    """
    buf = remote.get(off, off + 98304)
    i = buf.find(MAGIC)
    while i >= 0 and block_at(buf, i) is None:
        i = buf.find(MAGIC, i + 1)
    if i < 0:
        return None
    r = Reader(remote, off + len(buf), None)
    r.raw = buf[i:]
    for line in r.lines():
        return off + i, int(line.split(b"\t", 2)[1])
    return None


def fetch_chrom(chrom: str, positions: list[int], keep_cols: list[int], threads: int = 8) -> dict:
    """Rows of the dense file at the given positions, as raw VCF lines restricted to keep_cols sample columns."""
    url = SUBSET_URL.format(chrom=chrom)
    remote = Remote(url)
    tbi = CACHE / "tbi" / f"{chrom}.tbi"
    if not tbi.exists():
        tbi.parent.mkdir(parents=True, exist_ok=True)
        tbi.write_bytes(Remote(url + ".tbi").get(0, int(requests.head(url + ".tbi", timeout=60)
                                                         .headers["content-length"])))
    linear = read_tbi(tbi)[chrom]
    file_end = int(requests.head(url, timeout=60).headers["content-length"])
    by_window = defaultdict(list)
    for p in sorted(set(positions)):
        by_window[(p - 1) >> 14].append(p)
    out: dict[int, list[bytes]] = {}
    lock = threading.Lock()

    def window(w: int) -> None:
        start_v = linear[min(w, len(linear) - 1)]
        start = start_v >> 16
        later = [v >> 16 for v in linear[w + 1:w + 3] if v >> 16 > start]
        hi0 = later[-1] if later else file_end
        probes: list[tuple[int, int]] = []  # (block offset, first line POS), sorted by offset
        reader, last_pos, pending = None, -1, None
        for p in by_window[w]:
            lo, lo_skip = start, start_v & 0xFFFF
            hi = hi0
            for off, q in probes:
                if q < p and off > lo:
                    lo, lo_skip = off, None
                if q >= p and off < hi:
                    hi = off
            while hi - lo > 160 * 1024:
                mid = (lo + hi) // 2
                pr = probe(remote, mid)
                if pr is None or pr[0] >= hi:
                    hi = mid
                    continue
                bisect.insort(probes, pr)
                if pr[1] < p:
                    lo, lo_skip = pr[0], None
                else:
                    hi = mid
            # Reuse the open reader if it has not passed p and is no further back than lo.
            if reader is None or last_pos >= p or reader.position() < lo - 64 * 1024:
                reader, last_pos, pending = Reader(remote, lo, lo_skip), -1, None
                reader.it = reader.lines()
            rows = []
            while True:
                line = pending if pending is not None else next(reader.it, None)
                pending = None
                if line is None:
                    break
                q = int(line.split(b"\t", 2)[1])
                last_pos = q
                if q > p:
                    pending = line
                    last_pos = p
                    break
                if q == p:
                    f = line.split(b"\t")
                    rows.append(b"\t".join(f[:9] + [f[9 + i] for i in keep_cols]))
            with lock:
                out[p] = rows

    with ThreadPoolExecutor(threads) as ex:
        list(ex.map(window, sorted(by_window)))
    return {"rows": out, "bytes": remote.bytes, "requests": remote.requests}


def header(chrom: str) -> list[str]:
    """Sample names of the dense file (from its header)."""
    remote = Remote(SUBSET_URL.format(chrom=chrom))
    r = Reader(remote, 0, 0)
    for line in r.lines():
        if line.startswith(b"#CHROM"):
            return line.decode().split("\t")[9:]
    raise RuntimeError("no #CHROM line")


def slim_line(line: bytes) -> str:
    f = line.decode().split("\t")
    info = dict(kv.partition("=")[::2] for kv in f[7].split(";"))
    keep = [f"{k}={info[k]}" for k in INFO_KEEP if k in info]
    if "," in f[4]:
        keep.append("MULTI")
    fmt = f[8].split(":")
    idx = [fmt.index(x) if x in fmt else None for x in FMT_KEEP]
    gts = []
    for g in f[9:]:
        parts = g.split(":")
        gts.append(":".join(parts[i] if i is not None and i < len(parts) else "." for i in idx))
    return "\t".join(f[:6] + [".", ";".join(keep) or "."] + [":".join(FMT_KEEP)] + gts)


def fetch_one(chrom: str) -> None:
    """Fetch the rows at positions not fetched before; each run adds one part file (the cache is incremental)."""
    done_file = DENSE / f"{chrom}.positions.txt"
    done = set(map(int, done_file.read_text().split())) if done_file.exists() else set()
    positions = sorted({int(r["pos"]) for s in present_sets() for r in person_variants(s) if r["chrom"] == chrom} - done)
    if not positions:
        log(f"{chrom}: nothing new to fetch")
        return
    t0 = time.time()
    meta = samples()
    names = header(chrom)
    keep_cols = [i for i, s in enumerate(names) if meta[s]["kgp_release"] or meta[s]["subset"]]
    res = fetch_chrom(chrom, positions, keep_cols)
    part = len(list(DENSE.glob(f"{chrom}.part*.vcf.gz")))
    out = DENSE / f"{chrom}.part{part}.vcf.gz"
    tmp = DENSE / f"{chrom}.tmp.vcf"
    with tmp.open("w") as h:
        h.write("##fileformat=VCFv4.2\n")
        h.write(f"##contig=<ID={chrom}>\n")
        for k in INFO_KEEP:
            num = "A" if "_AC" in k or k.startswith("AC") else "1"
            h.write(f'##INFO=<ID={k},Number={num},Type=Integer,Description="{k} of the HGDP + 1KG subset file">\n')
        h.write('##INFO=<ID=MULTI,Number=0,Type=Flag,Description="Row had more than one ALT before splitting">\n')
        h.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        h.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        h.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        h.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        h.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                          + [names[i] for i in keep_cols]) + "\n")
        n = 0
        for p in sorted(res["rows"]):
            for line in res["rows"][p]:
                h.write(slim_line(line) + "\n")
                n += 1
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(out) + ".tmp", str(tmp)], check=True)
    Path(str(out) + ".tmp").rename(out)
    tmp.unlink()
    with done_file.open("a") as h:
        h.write("".join(f"{p}\n" for p in positions))
    (CACHE / f"counts_{chrom}.tsv").unlink(missing_ok=True)  # counts are redone from all parts
    hit = sum(1 for p in positions if res["rows"][p])
    log(f"{chrom}: {len(positions)} new positions, {hit} with a row, {n} rows; {res['requests']} requests, "
        f"{res['bytes'] / 1e6:.0f} MB in {time.time() - t0:.0f} s")


def step_fetch(chroms: list[str], procs: int = 14) -> None:
    DENSE.mkdir(parents=True, exist_ok=True)
    (CACHE / "logs").mkdir(exist_ok=True)
    running: list[tuple[str, subprocess.Popen]] = []
    queue = list(chroms)
    tries = Counter()
    while queue or running:
        while queue and len(running) < procs:
            c = queue.pop(0)
            tries[c] += 1
            running.append((c, subprocess.Popen([sys.executable, __file__, "fetch1", c],
                                                stdout=(CACHE / "logs" / f"{c}.log").open("a"),
                                                stderr=subprocess.STDOUT)))
        time.sleep(5)
        for c, p in list(running):
            if p.poll() is not None:
                running.remove((c, p))
                last = (CACHE / "logs" / f"{c}.log").read_text().splitlines()[-1:]
                if p.returncode == 0:
                    log(f"fetch: {last[0] if last else c}")
                elif tries[c] < 4:
                    log(f"fetch: {c} failed (exit {p.returncode}), retrying")
                    queue.append(c)
                else:
                    raise SystemExit(f"fetch: {c} failed {tries[c]} times")
    log("fetch: all chromosomes cached")


# ---------------------------------------------------------------------------------------------------------------
# Step 3: adj counts per row

def adj_call(g: str) -> tuple[int, int, int] | None:
    """(alt alleles, alleles, DP) of a call passing the adj filter, or None. g is GT:GQ:DP:AD after splitting.

    GQ at least 20, DP at least 10 (5 for a haploid call), and for a heterozygote AD of the alt allele over DP at
    least 0.2 (gnomAD's allele balance, get_adj_expr in gnomad_methods). The DP is returned so the caller can hold
    haploid calls to HAPLOID_DP.
    """
    gt, gq, dp, ad = (g.split(":") + [".", ".", "."])[:4]
    alleles = gt.replace("|", "/").split("/")
    if "." in alleles:
        return None
    try:
        gq_i, dp_i = int(gq), int(dp)
    except ValueError:
        return None
    haploid = len(alleles) == 1
    if gq_i < 20 or dp_i < (5 if haploid else 10):
        return None
    if not haploid and alleles[0] != alleles[1]:
        try:
            if int(ad.split(",")[1]) / dp_i < 0.2:
                return None
        except (ValueError, IndexError, ZeroDivisionError):
            return None
    return sum(a == "1" for a in alleles), len(alleles), dp_i


def count_one(chrom: str) -> None:
    meta = samples()
    fam_ids = {x for p in people().values() for x in p["family"]}
    out = CACHE / f"counts_{chrom}.tsv"
    tmp = CACHE / f"counts_{chrom}.tsv.tmp"
    parts = " ".join(str(p) for p in sorted(DENSE.glob(f"{chrom}.part*.vcf.gz")))
    q = subprocess.Popen(
        f"set -o pipefail; for p in {parts}; do bcftools norm -m- $p -Ou 2>/dev/null | bcftools query -H -f "
        f"'%CHROM\\t%POS\\t%REF\\t%ALT\\t%MULTI\\t" + "\\t".join(f"%INFO/{k}" for k in INFO_KEEP)
        + "[\\t%GT:%GQ:%DP:%AD]\\n' || exit 1; done", shell=True, executable="/bin/bash", stdout=subprocess.PIPE, text=True)
    cols = None
    kgp_cols = ["k_ac_" + g if i == 0 else "k_an_" + g for g in KGP_GROUPS for i in (0, 1)]
    with tmp.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(KEY + ["multi"] + INFO_KEEP + ["sub_ac", "sub_an", "sub_ac_h5", "sub_an_h5"] + kgp_cols
                   + ["fam"])
        for line in q.stdout:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                if cols is not None:
                    continue  # header of a later part; parts share the sample list
                names = subprocess.run(["bcftools", "query", "-l", parts.split()[0]], capture_output=True,
                                       text=True, check=True).stdout.split()
                cols = [(meta[s]["subset"], meta[s]["kgp_release"],
                         meta[s]["pop"] if meta[s]["pop"] in GROUPS else "other", s in fam_ids, s) for s in names]
                continue
            pos, ref, alt = trim(int(f[1]), f[2], f[3])
            info = [x if x != "." else "" for x in f[5:5 + len(INFO_KEEP)]]
            sub = [0, 0, 0, 0]
            kgp = {g: [0, 0] for g in KGP_GROUPS}
            fam = []
            for (in_sub, in_kgp, grp, in_fam, s), g in zip(cols, f[5 + len(INFO_KEEP):]):
                c = adj_call(g)
                if c is None:
                    continue
                ac, an, dp = c
                if in_sub:  # the protocol's DP 5 for haploid calls, kept only as a check
                    sub[2] += ac
                    sub[3] += an
                if an == 1 and dp < HAPLOID_DP:
                    continue
                if in_sub:
                    sub[0] += ac
                    sub[1] += an
                if in_kgp:
                    kgp[grp][0] += ac
                    kgp[grp][1] += an
                    if in_fam:
                        fam.append(f"{s}:{ac}:{an}")
            w.writerow([f[0], pos, ref, alt, 1 if f[4] == "1" else 0] + info + sub
                       + [x for g in KGP_GROUPS for x in kgp[g]] + [",".join(fam)])
    q.wait()
    assert q.returncode == 0, f"{chrom}: bcftools failed"
    tmp.rename(out)
    log(f"{chrom}: counts written")


def step_counts(chroms: list[str], procs: int = 12) -> None:
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(procs) as ex:
        list(ex.map(count_one, [c for c in chroms if not (CACHE / f"counts_{c}.tsv").exists()]))


def load_counts() -> dict[tuple, dict]:
    out = {}
    dup = 0
    for c in CHROMS:
        with (CACHE / f"counts_{c}.tsv").open() as h:
            for r in csv.DictReader(h, delimiter="\t"):
                k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                if k in out:
                    dup += 1
                    continue
                r["fam"] = {x.split(":")[0]: (int(x.split(":")[1]), int(x.split(":")[2]))
                            for x in r["fam"].split(",") if x}
                out[k] = r
    if dup:
        log(f"counts: {dup} duplicate keys after splitting (first row kept)")
    return out


# ---------------------------------------------------------------------------------------------------------------
# Step 4: leave-out frequencies

def leave_out(full: dict, sub_ac: dict, sub_an: dict) -> tuple[float, float, int]:
    """(AF, grpmax, number of groups where AC went negative) after subtracting per-group counts.

    sub_ac and sub_an are per group (GROUPS plus "other"); overall subtracts every group.
    """
    clamped = 0
    afs = {}
    for g in GROUPS:
        ac, an = full[f"ac_{g}"] - sub_ac.get(g, 0), full[f"an_{g}"] - sub_an.get(g, 0)
        if ac < 0:
            clamped += 1
            ac = 0
        afs[g] = ac / an if an > 0 else 0.0
    ac = full["ac"] - sum(sub_ac.values())
    an = full["an"] - sum(sub_an.values())
    if ac < 0:
        clamped += 1
        ac = 0
    return (ac / an if an > 0 else 0.0), max(afs.values()), clamped


def frequencies(key: tuple, members: list[str], full: dict | None, row: dict | None, meta: dict,
                stats: Counter) -> dict:
    """af_fam, grpmax_fam, af_no1kg, grpmax_no1kg and in_subset for one person variant."""
    chrom, pos = key[0], key[1]
    in_subset = int(row is not None)
    if full is None or not full["found"]:
        stats["not found in gnomAD"] += 1
        if row is not None:
            stats["in subset but not found in MyVariant"] += 1
        return {"af_fam": 0.0, "grpmax_fam": 0.0, "af_no1kg": 0.0, "grpmax_no1kg": 0.0, "in_subset": in_subset}
    grp = lambda s: meta[s]["pop"] if meta[s]["pop"] in GROUPS else "other"  # noqa: E731
    fam_ac, fam_an = Counter(), Counter()
    if row is not None:
        for s in members:
            ac, an = row["fam"].get(s, (0, 0))
            fam_ac[grp(s)] += ac
            fam_an[grp(s)] += an
        k_ac = {g: int(row[f"k_ac_{g}"]) for g in KGP_GROUPS}
        k_an = {g: int(row[f"k_an_{g}"]) for g in KGP_GROUPS}
    else:
        for s in members:
            fam_an[grp(s)] += 1 if haploid_region(chrom, pos, meta[s]["sex"]) else 2
        k_ac = {}
        k_an = {g: full[f"an_{g}"] * KGP_N[g] / N_GNOMAD[g] for g in GROUPS}
        k_an["other"] = full["an"] * 2500 / N_GNOMAD_ALL - sum(k_an.values())
    af_fam, gm_fam, c1 = leave_out(full, fam_ac, fam_an)
    af_no, gm_no, c2 = leave_out(full, k_ac, k_an)
    stats["family AC clamped at 0 (groups or overall)"] += c1
    stats["all-1KG AC clamped at 0 (groups or overall)"] += c2
    return {"af_fam": af_fam, "grpmax_fam": gm_fam, "af_no1kg": af_no, "grpmax_no1kg": gm_no, "in_subset": in_subset}


KGP_N: dict[str, int] = {}


def step_af() -> None:
    meta = samples()
    KGP_N.update(Counter(s["pop"] for s in meta.values() if s["kgp_release"]))
    full = load_full()
    counts = load_counts()
    ppl = people()
    all_stats = {}
    for set_name in present_sets():
        stats = Counter()
        rows = person_variants(set_name)
        out = DATA / set_name / "person_variant_af.tsv"
        with out.open("w", newline="") as h:
            w = csv.writer(h, delimiter="\t", lineterminator="\n")
            w.writerow(["sample"] + KEY + ["af_fam", "grpmax_fam", "af_no1kg", "grpmax_no1kg", "in_subset"])
            for r in rows:
                key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                f = frequencies(key, ppl[r["sample"]]["family"], full.get(key), counts.get(key), meta, stats)
                w.writerow([r["sample"], *key] + [f"{f[c]:.6g}" for c in ("af_fam", "grpmax_fam", "af_no1kg",
                                                                          "grpmax_no1kg")] + [f["in_subset"]])
                stats["rows"] += 1
                stats["in subset"] += f["in_subset"]
        log(f"af: {set_name}: {dict(stats)}")
        all_stats[set_name] = dict(stats)
    (CACHE / "af_stats.json").write_text(json.dumps(all_stats, indent=1))

# ---------------------------------------------------------------------------------------------------------------
# Step 5: checks

def sites_counts(keys: list[tuple]) -> dict[tuple, dict]:
    """AC and AN, overall and per group, from the gnomAD v3.1.2 genome sites VCF (the browser's source)."""
    out = {}

    def one(k: tuple) -> None:
        chrom, pos, ref, alt = k
        url = SITES_URL.format(chrom=chrom)
        tbi = CACHE / "tbi" / f"sites_{chrom}.tbi"
        if not tbi.exists():
            size = int(requests.head(url + ".tbi", timeout=60).headers["content-length"])
            tbi.write_bytes(Remote(url + ".tbi").get(0, size))
        fmt = "%REF\t%ALT\t%AC\t%AN\t" + "\t".join(f"%AC_{g}\t%AN_{g}" for g in GROUPS) + "\n"
        for _ in range(5):
            r = subprocess.run(["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}", "-f", fmt,
                                f"{url}##idx##{tbi}"], capture_output=True, text=True)
            if r.returncode == 0:
                break
            time.sleep(10)
        else:
            return
        for line in r.stdout.splitlines():
            f = line.split("\t")
            if (f[0], f[1]) == (ref, alt):
                vals = [int(x) if x != "." else 0 for x in f[2:]]
                out[k] = dict(zip(["ac", "an"] + [f"{x}_{g}" for g in GROUPS for x in ("ac", "an")], vals))

    (CACHE / "tbi").mkdir(parents=True, exist_ok=True)
    for c in sorted({k[0] for k in keys}):  # fetch each index once before the threads start
        one(next(k for k in keys if k[0] == c))
    with ThreadPoolExecutor(8) as ex:
        list(ex.map(one, [k for k in keys if k not in out]))
    return out


def pct(a: int, b: int) -> str:
    return f"{a:,} of {b:,} ({100 * a / b:.2f}%)" if b else "0 of 0"


def step_checks() -> None:
    meta = samples()
    KGP_N.update(Counter(s["pop"] for s in meta.values() if s["kgp_release"]))
    full = load_full()
    counts = load_counts()
    ppl = people()
    af_stats = json.loads((CACHE / "af_stats.json").read_text())
    md = ["# gnomAD frequency checks", "",
          "Made by `gnomad_af.py checks` (PROTOCOL.md section 1). Full counts: MyVariant.info `gnomad_genome` "
          "(gnomAD v3.1.2 genomes). Subset rows: `gnomad.genomes.v3.1.2.hgdp_tgp.chrN.vcf.bgz` at the positions of "
          "the 100 people's variants.", ""]

    # 1. The adj filter reproduces the subset file's own AC and AN.
    md += ["## 1. Adj counts against the subset file's own AC and AN", "",
           "The subset file's INFO AC and AN count the samples that are `high_quality` and not `related` in the "
           "sample metadata (3,380 of 4,151; 2,510 of them are 1000 Genomes samples, which is not the same set as the "
           "2,500 released ones, so the released samples' counts have no INFO field to compare with directly). "
           "Recounting those samples with the adj filter used here checks the filter and the haploid rule. Rows "
           "that were multiallelic before splitting are left out.", "",
           "| Rows | Checked | AC and AN both equal (adj as used here, haploid DP 10) | Same, with haploid DP 5 |",
           "|---|---|---|---|"]
    for label, sel in (("autosomes", lambda k: k[0] != "chrX"),
                       ("chrX outside the PARs", lambda k: k[0] == "chrX" and not in_par(k[0], k[1])),
                       ("chrX PARs", lambda k: k[0] == "chrX" and in_par(k[0], k[1])), ("all", lambda k: True)):
        rows = [r for k, r in counts.items() if sel(k) and r["multi"] == "0" and r["AC"] and r["AN"]]
        ok = sum(int(r["AC"]) == int(r["sub_ac"]) and int(r["AN"]) == int(r["sub_an"]) for r in rows)
        ok_dp = sum(int(r["AC"]) == int(r["sub_ac_h5"]) and int(r["AN"]) == int(r["sub_an_h5"]) for r in rows)
        md.append(f"| {label} | {len(rows):,} | {pct(ok, len(rows))} | {pct(ok_dp, len(rows))} |")
    bad = [(k, r) for k, r in counts.items() if r["multi"] == "0" and r["AC"]
           and (int(r["AC"]), int(r["AN"])) != (int(r["sub_ac"]), int(r["sub_an"]))][:5]
    if bad:
        md += ["", "First disagreements (key: INFO AC/AN, recounted AC/AN): " + "; ".join(
            f"{k[0]}:{k[1]} {k[2]}>{k[3]}: {r['AC']}/{r['AN']}, {r['sub_ac']}/{r['sub_an']}" for k, r in bad)]
    multi = sum(r["multi"] == "1" for r in counts.values())
    md += ["", "On chrX outside the PARs the subset file's own counts hold haploid calls to DP 10, not the "
           "protocol's 5, so the frequencies here use DP 10 for haploid calls (a deviation). The rules differ only "
           "for a haploid call with GQ 20 or more and DP 5 to 9, which changes a family member's AN by 1 allele."]
    md += ["", f"Rows that came from a multiallelic row: {multi:,} of {len(counts):,}. After splitting, a 1/2 "
           "heterozygote's allele balance uses the reference and that one alt allele only.", ""]

    # 2. gnomad_AC and gnomad_AN in the subset file against MyVariant's full counts.
    fields = ["ac", "an"] + [f"{x}_{g}" for g in GROUPS for x in ("ac", "an")]
    info_name = {"ac": "gnomad_AC", "an": "gnomad_AN", **{f"{x}_{g}": f"gnomad_{x.upper()}_{g}" for g in GROUPS
                                                         for x in ("ac", "an")}}
    both = [(k, r) for k, r in counts.items() if r["gnomad_AC"] and k in full and full[k]["found"]]
    ok = sum(all(int(r[info_name[f]] or 0) == full[k][f] for f in fields) for k, r in both)
    ok_all = sum(int(r["gnomad_AC"]) == full[k]["ac"] and int(r["gnomad_AN"]) == full[k]["an"] for k, r in both)
    no_info = sum(1 for r in counts.values() if not r["gnomad_AC"])
    md += ["## 2. Full counts: the subset file's gnomad_AC and gnomad_AN against MyVariant.info", "",
           f"- Rows with gnomad_ fields and a MyVariant.info record: {len(both):,}.",
           f"- Overall AC and AN equal: {pct(ok_all, len(both))}.",
           f"- Overall and all five groups' AC and AN equal: {pct(ok, len(both))}.",
           f"- Subset rows without gnomad_ fields (sites not in the full release): {no_info:,}.", ""]

    # 3. Released 1000 Genomes counts never exceed gnomAD's.
    over = Counter()
    for k, r in both:
        for g in GROUPS:
            if int(r[f"k_ac_{g}"]) > full[k][f"ac_{g}"]:
                over[f"AC {g}"] += 1
            if int(r[f"k_an_{g}"]) > full[k][f"an_{g}"]:
                over[f"AN {g}"] += 1
    md += ["## 3. Released 1000 Genomes adj counts within gnomAD's group counts", "",
           f"Rows where the 2,500 released samples' count exceeds the full gnomAD count in a group: "
           f"{dict(over) or 'none'} (of {len(both):,} rows).", ""]

    # 4. Group sizes used for the all-1KG fallback.
    md += ["## 4. Group sizes for the all-1000-Genomes fallback", "",
           "| Group | gnomAD v3.1.2 genomes (N) | Largest AN / 2 on autosomes in gnomad_full.tsv | Released 1000 Genomes |",
           "|---|---|---|---|"]
    for g in GROUPS:
        mx = max(v[f"an_{g}"] for k, v in full.items() if k[0] != "chrX") // 2
        md.append(f"| {g} | {N_GNOMAD[g]:,} | {mx:,} | {KGP_N[g]} |")
    md.append("")

    # 5. Matching and clamping.
    md += ["## 5. Matching, absent rows and clamping", ""]
    for set_name in present_sets():
        st = af_stats[set_name]
        md.append(f"- {set_name}: {st['rows']:,} person variants; in the subset file {pct(st.get('in subset', 0), st['rows'])}; "
                  f"not found in gnomAD (all frequencies 0) {st.get('not found in gnomAD', 0):,}; in the subset file "
                  f"but not found by MyVariant.info {st.get('in subset but not found in MyVariant', 0):,}; negative AC "
                  f"clamped at 0 (group or overall counts): family left out "
                  f"{st.get('family AC clamped at 0 (groups or overall)', 0):,}, all 1000 Genomes left out "
                  f"{st.get('all-1KG AC clamped at 0 (groups or overall)', 0):,}.")
    md.append("")

    # 6. Hand-worked examples.
    md += hand_examples(meta, full, counts, ppl)
    out = RESULTS / "gnomad_checks.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(md) + "\n")
    log(f"checks: wrote {out}")


def hand_examples(meta: dict, full: dict, counts: dict, ppl: dict) -> list[str]:
    """20 person variants (five chrX outside the PARs, one in a PAR if any), worked by hand, and a browser check."""
    rng = random.Random(1)
    pv = [(s, r) for s in present_sets() for r in person_variants(s)]
    keyed = [((r["chrom"], int(r["pos"]), r["ref"], r["alt"]), r["sample"]) for _, r in pv]
    keyed = [(k, s) for k, s in keyed if k in full and full[k]["found"]]
    xnon = [x for x in keyed if x[0][0] == "chrX" and not in_par(x[0][0], x[0][1])]
    xpar = [x for x in keyed if x[0][0] == "chrX" and in_par(x[0][0], x[0][1])]
    with_fam_row = [x for x in keyed if x[0] in counts and ppl[x[1]]["family"] and x[0][0] != "chrX"]
    no_row = [x for x in keyed if x[0] not in counts and ppl[x[1]]["family"] and x[0][0] != "chrX"]
    xnon_m = [x for x in xnon if ppl[x[1]]["sex"] == "male" and ppl[x[1]]["family"] and x[0] in counts]
    xnon_f = [x for x in xnon if ppl[x[1]]["sex"] == "female" and ppl[x[1]]["family"] and x[0] in counts]
    pick = rng.sample(xnon_m, min(3, len(xnon_m))) + rng.sample(xnon_f, min(2, len(xnon_f)))
    pick += rng.sample(xpar, min(1, len(xpar)))
    pick += rng.sample(no_row, min(2, len(no_row)))
    pick += rng.sample(with_fam_row, 20 - len(pick))
    sites = sites_counts(sorted({k for k, _ in pick}))
    lines = ["## 6. Hand-worked examples", "",
             f"Twenty person variants: {len([p for p in pick if p in xnon_m])} on chrX outside the PARs in XY people, "
             f"{len([p for p in pick if p in xnon_f])} in XX people, {len([p for p in pick if p in xpar])} in a PAR "
             f"({len(xpar)} person variants lie in a PAR), {len([p for p in pick if p in no_row])} with no subset row "
             "(fallback AN), the rest autosomal with a subset row. For each: the full counts from MyVariant.info and "
             "from the gnomAD v3.1.2 sites VCF (what the browser shows), each released family member's adj call "
             "(alt alleles / alleles), and the subtraction in the person's group and overall.", ""]
    agree = 0
    stats = Counter()
    for k, s in pick:
        f = full[k]
        st = sites.get(k)
        same = st is not None and all(st[x] == f[x] for x in st)
        agree += same
        p = ppl[s]
        row = counts.get(k)
        fam = []
        for m in p["family"]:
            g = meta[m]["pop"]
            if row is not None:
                ac, an = row["fam"].get(m, (0, 0))
                fam.append(f"{m} ({g}, {meta[m]['sex']}) {ac}/{an}")
            else:
                an = 1 if haploid_region(k[0], k[1], meta[m]["sex"]) else 2
                fam.append(f"{m} ({g}, {meta[m]['sex']}) no row, 0/{an}")
        res = frequencies(k, p["family"], f, row, meta, stats)
        grp = p["gnomad_pop"] if p["gnomad_pop"] in GROUPS else None
        lines.append(f"**{k[0]}:{k[1]} {k[2]}>{k[3]}**, {s} ({p['sex']}, gnomAD group {p['gnomad_pop'] or 'none'}, "
                     f"released {p['gnomad_release']}){', PAR' if in_par(k[0], k[1]) else ''}")
        lines.append(f"- Full: AC {f['ac']} / AN {f['an']}; " + ", ".join(
            f"{g} {f[f'ac_{g}']}/{f[f'an_{g}']}" for g in GROUPS) + f"; grpmax {f['grpmax']:.3g}. Sites VCF "
            + ("agrees." if same else f"differs or missing: {st}."))
        lines.append(f"- Family: {'; '.join(fam) or 'no released member'}.")
        fam_ac = sum(row["fam"].get(m, (0, 0))[0] for m in p["family"]) if row else 0
        fam_an = sum(row["fam"].get(m, (0, 0))[1] for m in p["family"]) if row else sum(
            1 if haploid_region(k[0], k[1], meta[m]["sex"]) else 2 for m in p["family"])
        work = f"overall ({f['ac']} - {fam_ac}) / ({f['an']} - {fam_an}) = {res['af_fam']:.4g}"
        if grp:
            g_ac = sum(row["fam"].get(m, (0, 0))[0] for m in p["family"] if meta[m]["pop"] == grp) if row else 0
            g_an = sum((row["fam"].get(m, (0, 0))[1] if row else
                        (1 if haploid_region(k[0], k[1], meta[m]["sex"]) else 2))
                       for m in p["family"] if meta[m]["pop"] == grp)
            a, n = max(f[f"ac_{grp}"] - g_ac, 0), f[f"an_{grp}"] - g_an
            work += f"; {grp} ({f[f'ac_{grp}']} - {g_ac}) / ({f[f'an_{grp}']} - {g_an}) = {a / n if n > 0 else 0:.4g}"
        lines.append(f"- Family left out: {work}; grpmax_fam {res['grpmax_fam']:.4g}. All 1000 Genomes left out: "
                     f"af {res['af_no1kg']:.4g}, grpmax {res['grpmax_no1kg']:.4g}.")
        lines.append("")
    lines.insert(4, f"Full counts agree with the sites VCF on {agree} of {len(pick)} (overall and five groups).")
    lines.insert(5, "")
    return lines


# ---------------------------------------------------------------------------------------------------------------
# Step 6: the issue #1 table and background sizes after the filter

def bins(afs: list[float]) -> list[str]:
    n = len(afs)
    b = [sum(a == 0 for a in afs), sum(0 < a < 1e-4 for a in afs), sum(1e-4 <= a < 1e-3 for a in afs),
         sum(a >= 1e-3 for a in afs)]
    return [f"{n:,}"] + [f"{100 * x / n:.0f}%" for x in b]


def step_overlap() -> None:
    full = load_full()
    with (DATA / "main" / "person_variant_af.tsv").open() as h:
        fam = {(r["sample"], r["chrom"], int(r["pos"]), r["ref"], r["alt"]): r for r in csv.DictReader(h, delimiter="\t")}
    with (V1_DATA / "cases.tsv").open() as h:
        cases = list(csv.DictReader(h, delimiter="\t"))
    key = lambda r: (r["chrom"], int(r["pos"]), r["ref"], r["alt"])  # noqa: E731
    spiked = {key(r) for r in cases if r["diagnostic"] == "1"}
    bg_rows = [r for r in cases if r["diagnostic"] == "0"]
    bg_keys = {key(r) for r in bg_rows}
    bg_pairs = {(r["sample"],) + key(r) for r in bg_rows}
    missing = [p for p in bg_pairs if p not in fam]
    assert not missing, f"{len(missing)} first-run background variants not in main/person_variant_af.tsv"
    faf = lambda k: full[k]["af"] if k in full else 0.0  # noqa: E731
    gmx = lambda k: full[k]["grpmax"] if k in full else 0.0  # noqa: E731
    md = ["# gnomAD overlap of the first run's 200 SNV cases, with the new frequencies", "",
          "Made by `gnomad_af.py overlap`. Issue #1's table, redone. Bins use the overall gnomAD v3.1.2 genome AF; "
          "\"not in gnomAD\" means AF 0 (not found, or no allele left after leaving people out). The spiked-in "
          "(diagnostic) variants are not carried by the person, so they keep the full frequency. Full-frequency rows "
          "count distinct variants, as in the issue; leave-out rows count person variants (the same variant in two "
          "people can have different frequencies).", "",
          "| Variants | Frequency | n | Not in gnomAD | Below 0.01% | 0.01 to 0.1% | 0.1% or more |",
          "|---|---|---|---|---|---|---|",
          "| 200 cases, spiked-in | full | " + " | ".join(bins([faf(k) for k in spiked])) + " |",
          "| 200 cases, background | full | " + " | ".join(bins([faf(k) for k in bg_keys])) + " |",
          "| 200 cases, background, per person variant | full | " + " | ".join(
              bins([faf(p[1:]) for p in bg_pairs])) + " |",
          "| 200 cases, background | family left out | " + " | ".join(
              bins([float(fam[p]["af_fam"]) for p in bg_pairs])) + " |",
          "| 200 cases, background | all 1000 Genomes left out | " + " | ".join(
              bins([float(fam[p]["af_no1kg"]) for p in bg_pairs])) + " |", ""]

    per_case = defaultdict(list)
    for r in bg_rows:
        per_case[r["case"]].append((r["sample"],) + key(r))
    md += ["## Background variants per case left after the grpmax filter", "",
           "Median (and range) over the 200 cases of the background SNVs kept, by the frequency used for grpmax.", "",
           "| Frequency | Before the filter | grpmax below 0.1% | grpmax below 1% |", "|---|---|---|---|"]
    for label, f in (("full", lambda p: gmx(p[1:])), ("family left out", lambda p: float(fam[p]["grpmax_fam"])),
                     ("all 1000 Genomes left out", lambda p: float(fam[p]["grpmax_no1kg"]))):
        cells = []
        for cut in (None, 1e-3, 1e-2):
            n = [sum(1 for p in v if cut is None or f(p) < cut) for v in per_case.values()]
            cells.append(f"{statistics.median(n):g} ({min(n)} to {max(n)})")
        md.append(f"| {label} | " + " | ".join(cells) + " |")
    md += ["", f"All {len(bg_pairs):,} background person variants of the 200 cases join to "
           "`data/main/person_variant_af.tsv`."]
    out = RESULTS / "gnomad_overlap.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(md) + "\n")
    log(f"overlap: wrote {out}")
    print("\n".join(md))


# ---------------------------------------------------------------------------------------------------------------
# Command line

def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    chroms = sys.argv[2:] or CHROMS
    if cmd == "fetch1":
        fetch_one(chroms[0])
        return
    steps = {"full": step_full, "fetch": lambda: step_fetch(chroms), "counts": lambda: step_counts(chroms),
             "af": step_af, "checks": step_checks, "overlap": step_overlap}
    for name in (list(steps) if cmd == "all" else [cmd]):
        log(f"== {name}")
        steps[name]()


if __name__ == "__main__":
    main()
