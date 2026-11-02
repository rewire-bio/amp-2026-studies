"""Look up precomputed scores for every case variant, the way a clinical pipeline would.

Implements PROTOCOL.md sections 4 (CADD consequence rule), 6 (scores) and 10 (EVEE join).
Derived from rare-disease/lookup_scores.py at amp-2026-v4.

  uv run -q --python 3.12 --with-requirements requirements.txt python lookup_scores.py --study snv
  uv run ... python lookup_scores.py --study indel                       # adds spliceai_run
  uv run ... python lookup_scores.py --study replication cadd avi        # some sources only
  uv run ... python lookup_scores.py --variants pool.tsv --out-dir data/pool cadd   # any variant list

Reads data/<study>/case_variants.tsv (or --variants) and writes scores_<source>.tsv next to it (or in --out-dir),
one row per variant found, keyed by chrom, pos, ref, alt (chr-prefixed), sorted by the key. A variant with no row
was looked up and not found.

Sources (GRCh38):
  cadd          CADD v1.7 inclAnno by remote tabix (common.cadd_lookup): cadd_phred; cadd_consequence, the most
                severe consequence over MANE Select rows (all rows if none is MANE); cadd_consequence_any, the most
                severe over all rows; gene; mane (1 if a MANE row exists, else 0); spliceai_max (largest of the four
                SpliceAI deltas). Indels come from CADD's gnomAD v4.0 indel table only.
  revel         REVEL v1.3 zip, max over transcripts. SNVs only.
  am            AlphaMissense hg38 table (../rare-disease/data/am), max over transcripts. SNVs only.
  avi           AlphaGenome Variant Impact all-SNV table (../rare-disease/data/avi, a link to an external SSD) by
                local tabix: avi (raw score) and avi_phred. SNVs only.
  gpnmsa        GPN-MSA all-SNV table (Hugging Face songlab/gpn-msa-hg38-scores) by remote tabix; score is -LLR.
  gpnstar       GPN-Star vertebrate (v100) and mammal (m447) models (songlab/gpn-star-scores Parquet) by remote
                DuckDB; -llr_calibrated. Writes scores_gpnstar_v.tsv and scores_gpnstar_m.tsv.
  spliceai_run  SpliceAI 1.3.1 deltas from data/spliceai/indels.out.vcf (run_spliceai.sh): max of the four.
  evee          EVEE (Zenodo 19701997): evee (pathogenicity), evee_label, evee_variation_id. EVEE positions are
                0-based, so a variant joins at EVEE pos equal to VCF POS minus 1. The slim table
                data/evee/evee_slim.parquet is built on first use from the five Zenodo shards (ten columns, byte
                ranges), and results/evee_source.md records the record version, file sizes, checksums and date.

Reuse, so that nothing already fetched is fetched again. Two sources of earlier rows:
  - The first run (REVEL, AlphaMissense, AVI, GPN-MSA and GPN-Star, which do not depend on annotation; not CADD,
    whose consequence rule changed). Every first-run directory (../rare-disease/data/{.,indel,replication,cancer})
    that has scores_<source>.tsv counts all keys of its case_variants.tsv as looked up: a key with a row is copied, a
    key without one was not found.
  - This run (every source above and CADD): the rows in data/snv and data/indel scores_<source>.tsv, and for CADD
    also data/pool/scores_cadd.tsv (pool.py, where mane is True or False; it is written here as 1 or 0). Only keys
    with a row count as looked up, because this run's case_variants.tsv may have been rebuilt since the scores were
    written; a key looked up before and not found is looked up again.
Only the remaining keys are looked up. Where a key has rows in more than one place, this run's row is kept, and the
number of keys whose values differ is printed. spliceai_run reuse happens in run_spliceai.sh. EVEE is a local join
and is not reused. --no-reuse looks everything up again (for tests).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

import common
from common import DATA, RESULTS, V1_DATA, is_snv

AM = V1_DATA / "am/AlphaMissense_hg38.tsv.gz"
AVI = V1_DATA / "avi/alphagenome_variant_impact_score_snvs.tsv.gz"
GPN_MSA = "https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz"
GPN_STAR = "hf://datasets/songlab/gpn-star-scores/data/gpn-star-hg38-{model}-200m/llr/llr_{chrom}.parquet"
SPLICEAI_VCF = DATA / "spliceai/indels.out.vcf"
EVEE_RECORD = "19701997"
EVEE_SHARD = "https://zenodo.org/records/" + EVEE_RECORD + "/files/clean_shard_{i}.parquet?download=1"
EVEE_COLUMNS = ["chrom", "pos", "ref", "alt", "pathogenicity", "label", "significance", "stars", "consequence",
                "variation_id"]
EVEE_SLIM = DATA / "evee/evee_slim.parquet"
EVEE_SOURCE = RESULTS / "evee_source.md"
V1_DIRS = [V1_DATA, V1_DATA / "indel", V1_DATA / "replication", V1_DATA / "cancer"]
# This run's earlier scores (rows only count as looked up). CADD also reuses pool.py's lookup of the ClinVar pool.
V2_DIRS = [DATA / "snv", DATA / "indel"]
V2_EXTRA = {"cadd": [DATA / "pool" / "scores_cadd.tsv"]}
# Local copies of tabix indexes the first run already downloaded (read only), so remote tabix need not fetch them.
V1_INDEXES = [V1_DATA / "whole_genome_SNVs_inclAnno.tsv.gz.tbi", V1_DATA / "gnomad.genomes.r4.0.indel_inclAnno.tsv.gz.tbi",
              V1_DATA / "gpn_msa.tsv.bgz.tbi"]

Key = tuple  # (chrom, pos, ref, alt)


def read_keys(path: Path) -> set[Key]:
    with path.open() as h:
        return {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in csv.DictReader(h, delimiter="\t")}


def first_run_rows(name: str, header: list[str]) -> tuple[set[Key], dict[Key, list[str]]]:
    """Keys the first run looked up for scores_<name>.tsv, and the rows it found."""
    looked, rows, clashes = set(), {}, 0
    for d in V1_DIRS:
        scores = d / f"scores_{name}.tsv"
        if not scores.exists():
            continue
        keys = read_keys(d / "case_variants.tsv")
        looked |= keys
        with scores.open() as h:
            for r in csv.DictReader(h, delimiter="\t"):
                k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                v = [r[c] for c in header]
                if k in rows and rows[k] != v:
                    clashes += 1
                rows[k] = v
    if clashes:
        print(f"  warning: {clashes} first-run {name} rows differ between directories; the last one read is kept")
    return looked, rows


def norm_value(col: str, v: str) -> str:
    """pool.py writes mane as True or False and spliceai_max as 0 where the study files write 1 or 0 and 0.0."""
    if col == "mane":
        return {"True": "1", "False": "0"}.get(v, v)
    if col == "spliceai_max" and v:
        return str(float(v))
    return v


def this_run_rows(name: str, header: list[str]) -> dict[Key, list[str]]:
    """Rows this run already found for scores_<name>.tsv (data/snv, data/indel, and the pool for CADD)."""
    rows, clashes = {}, 0
    files = [d / f"scores_{name}.tsv" for d in V2_DIRS] + V2_EXTRA.get(name, [])
    for scores in files:
        if not scores.exists():
            continue
        with scores.open() as h:
            for r in csv.DictReader(h, delimiter="\t"):
                k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                v = [norm_value(c, r[c]) for c in header]
                if k in rows and rows[k] != v:
                    clashes += 1
                rows.setdefault(k, v)
    if clashes:
        print(f"  warning: {clashes} {name} rows of this run differ between files; the first one read is kept")
    return rows


class Run:
    def __init__(self, wanted: set[Key], out: Path, reuse: bool, work: Path):
        self.wanted, self.out, self.reuse, self.work = wanted, out, reuse, work

    def split(self, name: str, header: list[str], first_run: bool = True,
              snv_only: bool = True) -> tuple[dict[Key, list[str]], set[Key]]:
        """Rows copied from earlier lookups (this run's, then the first run's), and the keys still to look up.

        snv_only leaves indels out of the keys to look up, for sources that score SNVs only.
        """
        if not self.reuse:
            return {}, snvs(self.wanted) if snv_only else set(self.wanted)
        looked, rows = first_run_rows(name, header) if first_run else (set(), {})
        own = this_run_rows(name, header)
        both = self.wanted & set(own) & set(rows)
        differ = sum(1 for k in both if own[k] != rows[k])
        copied = {k: rows[k] for k in self.wanted & looked if k in rows}
        copied |= {k: own[k] for k in self.wanted & set(own)}
        todo = self.wanted - looked - set(own)
        if snv_only:
            todo = snvs(todo)
        print(f"{name}: {len(self.wanted & set(own))} found by this run, {len(self.wanted & looked)} looked up by the "
              f"first run ({sum(1 for k in self.wanted & looked if k in rows)} found)"
              + (f", {differ} of {len(both)} in both differ (this run's kept)" if both else "")
              + f"; {len(todo)} to look up")
        return copied, todo

    def write(self, name: str, header: list[str], rows: dict[Key, list]) -> None:
        out = self.out / f"scores_{name}.tsv"
        with out.open("w", newline="") as h:
            w = csv.writer(h, delimiter="\t", lineterminator="\n")
            w.writerow(["chrom", "pos", "ref", "alt", *header])
            for k in sorted(rows):
                w.writerow([*k, *rows[k]])
        print(f"{name}: {len(rows)} of {len(self.wanted)} variants scored, wrote {out}", flush=True)


def snvs(keys: set[Key]) -> set[Key]:
    return {k for k in keys if is_snv(k[2], k[3])}


def tabix(args: list[str], cwd: Path, what: str, tries: int = 8) -> str:
    for attempt in range(tries):
        r = subprocess.run(["tabix", *args], capture_output=True, text=True, cwd=cwd)
        if r.returncode == 0:
            return r.stdout
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{what} tabix failed: {r.stderr[-300:]}")


# ---------------------------------------------------------------------------------------------------------------

def cadd(run: Run) -> None:
    header = ["cadd_phred", "cadd_consequence", "cadd_consequence_any", "gene", "mane", "spliceai_max"]
    rows, todo = run.split("cadd", header, first_run=False, snv_only=False)
    found = common.cadd_lookup(todo, run.work) if todo else {}
    rows |= {k: [d["cadd_phred"], d["cadd_consequence"], d["cadd_consequence_any"], d["gene"], int(d["mane"]),
                 "" if d["spliceai_max"] is None else d["spliceai_max"]] for k, d in found.items()}
    no_mane = sum(1 for v in rows.values() if str(v[4]) == "0")
    differs = sum(1 for v in rows.values() if str(v[4]) == "1" and v[1] != v[2])
    print(f"  CADD: {len(found)} of {len(todo)} looked up found; {len(rows)} rows in all, {no_mane} with no MANE row "
          f"(all-transcript consequence used), {differs} where the MANE consequence differs from the all-transcript one")
    run.write("cadd", header, rows)


def revel(run: Run) -> None:
    rows, todo = run.split("revel", ["revel"])
    todo = snvs(todo)
    if todo:
        new: dict[Key, float] = {}
        with zipfile.ZipFile(common.REVEL_ZIP) as z, z.open("revel_with_transcript_ids") as raw:
            for line in io.TextIOWrapper(raw):
                f = line.split(",")
                if f[2] == "." or not f[2].isdigit():
                    continue
                key = ("chr" + f[0], int(f[2]), f[3], f[4])
                if key in todo:
                    s = float(f[7])
                    new[key] = max(s, new.get(key, s))
        rows |= {k: [v] for k, v in new.items()}
    run.write("revel", ["revel"], rows)


def alphamissense(run: Run) -> None:
    rows, todo = run.split("am", ["alphamissense"])
    todo = snvs(todo)
    if todo:
        new: dict[Key, float] = {}
        with gzip.open(AM, "rt") as h:
            for line in h:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                key = (f[0], int(f[1]), f[2], f[3])
                if key in todo:
                    s = float(f[8])
                    new[key] = max(s, new.get(key, s))
        rows |= {k: [v] for k, v in new.items()}
    run.write("am", ["alphamissense"], rows)


def avi(run: Run) -> None:
    rows, todo = run.split("avi", ["avi", "avi_phred"])
    todo = snvs(todo)
    if todo:
        if not (AVI.exists() and Path(str(AVI) + ".tbi").exists()):
            raise SystemExit(f"AVI table not found at {AVI}. {V1_DATA / 'avi'} links to the external SSD "
                             f"(/Volumes/Extreme SSD/avi); mount it and run again. {len(todo)} SNVs need AVI.")
        regions = run.work / "avi_regions.tsv"
        regions.write_text("".join(f"{c}\t{p}\t{p}\n" for c, p in sorted({(c, p) for c, p, _, _ in todo})))
        found = 0
        for line in tabix(["-R", str(regions), str(AVI)], run.work, "AVI").splitlines():
            f = line.split("\t")
            key = (f[0], int(f[1]), f[2], f[3])
            if key in todo:
                rows[key] = [f[4], f[5]]
                found += 1
        print(f"  AVI: {found} of {len(todo)} SNVs found")
    run.write("avi", ["avi", "avi_phred"], rows)


def gpn_msa(run: Run) -> None:
    """tabix follows Hugging Face's redirect for the data but not for the index, so the index is fetched first."""
    rows, todo = run.split("gpnmsa", ["gpn_msa"])
    todo = snvs(todo)
    if todo:
        idx = run.work / "gpn_msa.tsv.bgz.tbi"
        if not idx.exists():
            subprocess.run(["curl", "-sSfL", "--retry", "5", "-o", str(idx), GPN_MSA + ".tbi"], check=True)

        def one_chrom(chrom: str) -> dict:
            regions = run.work / f"gpn_msa_regions_{chrom}.tsv"
            regions.write_text("".join(f"{chrom[3:]}\t{p}\t{p}\n"
                                       for p in sorted({p for c, p, _, _ in todo if c == chrom})))
            found = {}
            for line in tabix(["-R", str(regions), f"{GPN_MSA}##idx##{idx}"], run.work, f"GPN-MSA {chrom}").splitlines():
                f = line.split("\t")
                key = ("chr" + f[0], int(f[1]), f[2], f[3])
                if key in todo:
                    found[key] = [-float(f[4])]
            return found

        with ThreadPoolExecutor(4) as pool:
            for found in pool.map(one_chrom, sorted({c for c, *_ in todo})):
                rows |= found
    run.write("gpnmsa", ["gpn_msa"], rows)


def gpn_star(run: Run) -> None:
    for model, name in (("v100", "gpnstar_v"), ("m447", "gpnstar_m")):
        rows, todo = run.split(name, [name])
        todo = snvs(todo)
        if todo:
            import duckdb
            con = duckdb.connect()
            con.sql("INSTALL httpfs; LOAD httpfs;")
            for chrom in sorted({c for c, *_ in todo}):
                positions = sorted({p for c, p, _, _ in todo if c == chrom})
                q = (f"SELECT pos, ref, alt, llr_calibrated FROM read_parquet('{GPN_STAR.format(model=model, chrom=chrom)}') "
                     f"WHERE pos IN ({','.join(map(str, positions))})")
                for attempt in range(6):
                    try:
                        result = con.sql(q).fetchall()
                        break
                    except Exception as e:  # network errors from httpfs
                        print(f"  {name} {chrom}: retry after {type(e).__name__}", flush=True)
                        time.sleep(10 * (attempt + 1))
                else:
                    raise RuntimeError(f"GPN-Star {model} {chrom} failed")
                for pos, ref, alt, llr in result:
                    key = (chrom, pos, ref, alt)
                    if key in todo and llr is not None:
                        rows[key] = [-llr]
        run.write(name, [name], rows)


def spliceai_run(run: Run, vcf: Path = SPLICEAI_VCF) -> None:
    if not vcf.exists():
        raise SystemExit(f"{vcf} not found: run run_spliceai.sh first")
    rows, in_vcf = {}, set()
    with vcf.open() as h:
        for line in h:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]), f[3], f[4])
            in_vcf.add(key)
            if key not in run.wanted or not f[7].startswith("SpliceAI="):
                continue
            # One entry per gene: ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL; keep the largest delta.
            deltas = [float(x) for e in f[7][len("SpliceAI="):].split(",") for x in e.split("|")[2:6] if x != "."]
            if deltas:
                rows[key] = [max(deltas)]
    # run_spliceai.sh writes one row per input indel (INFO '.' when SpliceAI gives no score), so a wanted indel with
    # no row means the VCF is stale, not that the indel is unscored.
    stale = {k for k in run.wanted if not is_snv(k[2], k[3])} - in_vcf
    if stale:
        raise SystemExit(f"{len(stale)} indels have no row in {vcf}: rerun run_spliceai.sh")
    run.write("spliceai_run", ["spliceai_run"], rows)


# ---------------------------------------------------------------------------------------------------------------
# EVEE

def evee_build() -> None:
    """Read ten columns of the five Zenodo shards by byte range and write data/evee/evee_slim.parquet."""
    import fsspec
    import pyarrow as pa
    import pyarrow.parquet as pq
    EVEE_SLIM.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{EVEE_RECORD}", timeout=60) as r:
        record = json.load(r)

    def shard(i: int):
        # Each shard's slim part is kept until the whole table is written, so a failed build resumes.
        part = EVEE_SLIM.parent / f"evee_slim_part{i}.parquet"
        if part.exists():
            return pq.read_table(part)
        for attempt in range(6):
            try:
                t = time.time()
                fs = fsspec.filesystem("http")
                with fs.open(EVEE_SHARD.format(i=i), "rb", block_size=8 << 20) as f:
                    table = pq.ParquetFile(f).read(columns=EVEE_COLUMNS)
                pq.write_table(table, part)
                print(f"  EVEE shard {i}: {table.num_rows} rows, {time.time() - t:.0f} s", flush=True)
                return table
            except Exception as e:  # network errors
                print(f"  EVEE shard {i}: retry after {type(e).__name__}: {e}", flush=True)
                time.sleep(30 * (attempt + 1))
        raise RuntimeError(f"EVEE shard {i} failed")

    with ThreadPoolExecutor(3) as pool:  # Zenodo slows down with more parallel readers
        tables = list(pool.map(shard, range(5)))
    table = pa.concat_tables(tables)
    import duckdb
    con = duckdb.connect()
    con.register("slim", table)
    distinct = con.execute("SELECT count(DISTINCT (chrom, pos, ref, alt)) FROM slim").fetchall()[0][0]
    if distinct != table.num_rows:
        raise ValueError(f"EVEE keys are not unique: {table.num_rows} rows, {distinct} distinct keys")
    tmp = EVEE_SLIM.with_suffix(".tmp")
    pq.write_table(table, tmp)
    tmp.rename(EVEE_SLIM)
    for i in range(5):
        (EVEE_SLIM.parent / f"evee_slim_part{i}.parquet").unlink(missing_ok=True)
    shards =sorted((f for f in record["files"] if f["key"].startswith("clean_shard_")), key=lambda f: f["key"])
    meta = record["metadata"]
    lines = [
        "# EVEE source",
        "",
        "Written by `lookup_scores.py` when it built `data/evee/evee_slim.parquet` (PROTOCOL.md section 10).",
        "",
        f"- Record: https://zenodo.org/records/{EVEE_RECORD} ({meta.get('title', '')}), licence "
        f"{(meta.get('license') or {}).get('id', '')}",
        f"- Version: {meta.get('version') or 'not stated'}; revision {record.get('revision')}; published "
        f"{meta.get('publication_date')}; last updated {record.get('updated')}",
        f"- Downloaded: {date.today().isoformat()}, columns {', '.join(EVEE_COLUMNS)} only, read by byte range",
        f"- Rows: {table.num_rows} ({', '.join(str(t.num_rows) for t in tables)} per shard), all keys distinct",
        f"- `evee_slim.parquet` SHA-256: {hashlib.sha256(EVEE_SLIM.read_bytes()).hexdigest()}",
        "",
        "Only the ten columns were read, so the checksums below are Zenodo's for the full shards, not checked here.",
        "",
        "| File | Size (bytes) | Checksum (Zenodo) |",
        "|---|---|---|",
        *[f"| {f['key']} | {f['size']} | {f['checksum']} |" for f in shards],
        "",
    ]
    EVEE_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    EVEE_SOURCE.write_text("\n".join(lines))
    print(f"  wrote {EVEE_SLIM} ({table.num_rows} rows) and {EVEE_SOURCE}")


def evee(run: Run) -> None:
    import duckdb
    if not EVEE_SLIM.exists():
        evee_build()
    con = duckdb.connect()
    con.execute("CREATE TABLE w (chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR)")
    con.executemany("INSERT INTO w VALUES (?, ?, ?, ?)", sorted(run.wanted))
    # EVEE positions are 0-based: VCF POS minus 1. evee_build checks that keys are unique, so the join gives at most
    # one row per variant.
    result = con.execute(f"""
        SELECT w.chrom, w.pos, w.ref, w.alt, e.pathogenicity, e.label, e.variation_id
        FROM w JOIN read_parquet('{EVEE_SLIM}') e
          ON e.chrom = w.chrom AND e.pos = w.pos - 1 AND e.ref = w.ref AND e.alt = w.alt""").fetchall()
    rows: dict[Key, list] = {}
    for c, p, r, a, path, label, vid in result:
        if (c, p, r, a) in rows:
            raise ValueError(f"more than one EVEE row for {c}:{p} {r}>{a}")
        rows[(c, p, r, a)] = ["" if path is None else path, label or "", "" if vid is None else vid]
    run.write("evee", ["evee", "evee_label", "evee_variation_id"], rows)


# ---------------------------------------------------------------------------------------------------------------

SOURCES = {"cadd": cadd, "revel": revel, "am": alphamissense, "avi": avi, "gpnmsa": gpn_msa, "gpnstar": gpn_star,
           "spliceai_run": spliceai_run, "evee": evee}
DEFAULT = ["cadd", "revel", "am", "avi", "gpnmsa", "gpnstar", "evee"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("sources", nargs="*", help=f"any of {', '.join(SOURCES)}; default: all (spliceai_run for indel only)")
    p.add_argument("--study", choices=list(common.STUDIES), help="read data/<study>/case_variants.tsv")
    p.add_argument("--variants", type=Path, help="any table with chrom, pos, ref, alt (with --out-dir)")
    p.add_argument("--out-dir", type=Path)
    p.add_argument("--no-reuse", action="store_true", help="look everything up again (tests)")
    a = p.parse_args()
    if unknown := [s for s in a.sources if s not in SOURCES]:
        p.error(f"unknown sources: {', '.join(unknown)}")
    if bool(a.study) == bool(a.variants):
        p.error("give --study or --variants (with --out-dir)")
    if a.study:
        out, variants = common.study_dir(a.study), common.study_dir(a.study) / "case_variants.tsv"
    else:
        if not a.out_dir:
            p.error("--variants needs --out-dir")
        out, variants = a.out_dir, a.variants
    if not variants.exists():
        raise SystemExit(f"{variants} not found (build_cases.py writes the case variant lists)")
    out.mkdir(parents=True, exist_ok=True)
    wanted = read_keys(variants)
    print(f"{len(wanted)} variants ({len(snvs(wanted))} SNVs) from {variants}", flush=True)
    sources = a.sources or DEFAULT + (["spliceai_run"] if a.study == "indel" else [])
    with tempfile.TemporaryDirectory(prefix="lookup_scores_") as tmp:
        work = Path(tmp)
        for idx in V1_INDEXES:
            if idx.exists():
                shutil.copy(idx, work / idx.name)
        run = Run(wanted, out, not a.no_reuse, work)
        for name in sources:
            t = time.time()
            SOURCES[name](run)
            print(f"  ({name}: {time.time() - t:.0f} s)", flush=True)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
