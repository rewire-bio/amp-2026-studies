"""Look up precomputed scores for every case variant, the way a clinical pipeline would.

  python rare-disease/lookup_scores.py                       # all sources, SNV study (data/)
  python rare-disease/lookup_scores.py cadd                  # one source
  python rare-disease/lookup_scores.py --study indel         # data/indel/
  python rare-disease/lookup_scores.py --study replication --reuse   # copy rows for variants already in data/

Sources (GRCh38):
  cadd      CADD v1.7 inclAnno by remote tabix: PHRED, consequence, SpliceAI delta scores (max of four). Indels
            from CADD's precomputed gnomAD v4.0 indel table; indels not in it stay unscored, as in a lookup
  revel     REVEL v1.3 zip from evo2-cloud-operations, max over transcripts
  am        AlphaMissense hg38 table (data/am), max over transcripts
  gpnmsa    GPN-MSA all-SNV table (Hugging Face songlab/gpn-msa-hg38-scores) by remote tabix; score is -LLR
  spliceai_run  SpliceAI 1.3.1 run on the study's indels by run_spliceai.sh (max of the four delta scores)
  gpnstar   GPN-Star vertebrate and mammal models (songlab/gpn-star-scores Parquet) by remote DuckDB; -llr_calibrated
  avi       AlphaGenome Variant Impact score, static all-SNV table (avi_scores_snvs_tabix.zip, 88.5 GB, unpacked
            into data/avi/) by local tabix; raw_score. SNVs only

Outputs data/scores_<source>.tsv keyed by chrom, pos, ref, alt (chr-prefixed). Score tables are not committed.
"""
import argparse
import csv
import gzip
import io
import subprocess
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
CADD = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz"
CADD_INDEL = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel_inclAnno.tsv.gz"
OUT = DATA  # set by --study
REUSE: set = set()  # variants already looked up for the SNV study (--reuse)
REVEL_ZIP = HERE.parents[1] / "evo2-cloud-operations/data/raw/current-predictors/revel-v1.3_all_chromosomes.zip"
AM = DATA / "am/AlphaMissense_hg38.tsv.gz"
GPN_MSA = "https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz"
AVI = DATA / "avi/alphagenome_variant_impact_score_snvs.tsv.gz"
GPN_STAR = "hf://datasets/songlab/gpn-star-scores/data/gpn-star-hg38-{model}-200m/llr/llr_{chrom}.parquet"
SEVERITY = {c: i for i, c in enumerate(["STOP_GAINED", "STOP_LOST", "CANONICAL_SPLICE", "NON_SYNONYMOUS", "SPLICE_SITE",
                                        "SYNONYMOUS", "5PRIME_UTR", "3PRIME_UTR", "NONCODING_CHANGE", "INTRONIC",
                                        "UPSTREAM", "DOWNSTREAM", "REGULATORY", "INTERGENIC"])}


def variants(study: Path) -> set[tuple[str, int, str, str]]:
    with (study / "case_variants.tsv").open() as h:
        return {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in csv.DictReader(h, delimiter="\t")}


def reused(name: str, header: list[str]) -> dict:
    """Rows of data/scores_<name>.tsv for the variants in REUSE; a variant with no row there was not found."""
    if not REUSE:
        return {}
    with (DATA / f"scores_{name}.tsv").open() as h:
        return {k: [r[c] for c in header] for r in csv.DictReader(h, delimiter="\t")
                if (k := (r["chrom"], int(r["pos"]), r["ref"], r["alt"])) in REUSE}


def write(name: str, header: list[str], rows: dict) -> None:
    rows = reused(name, header) | rows
    out = OUT / f"scores_{name}.tsv"
    with out.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", *header])
        for k in sorted(rows, key=lambda k: (k[0], k[1])):
            w.writerow([*k, *rows[k]])
    print(f"{name}: {len(rows)} variants scored, wrote {out.name}")


def cadd(wanted: set) -> None:
    rows = {}
    snvs = {k for k in wanted if len(k[2]) == len(k[3]) == 1}
    for url, subset, tag in ((CADD, snvs, "snv"), (CADD_INDEL, wanted - snvs, "indel")):
        if subset:
            rows.update(cadd_file(url, subset, tag))
    write("cadd", ["cadd_phred", "cadd_consequence", "gene", "spliceai_max"],
          {k: [v[0], v[1], v[2], "" if v[3] is None else v[3]] for k, v in rows.items()})


def cadd_file(url: str, wanted: set, tag: str) -> dict:
    header = subprocess.check_output(["tabix", "-H", url], text=True, cwd=DATA).splitlines()[-1].lstrip("#")
    col = {c: i for i, c in enumerate(header.split("\t"))}
    splice = [col[c] for c in ("SpliceAI-acc-gain", "SpliceAI-acc-loss", "SpliceAI-don-gain", "SpliceAI-don-loss")]

    def one_chrom(chrom: str) -> dict:
        positions = sorted({p for c, p, _, _ in wanted if c == chrom})
        regions = DATA / f"kgp/cadd_regions_{tag}_{chrom}.tsv"
        regions.write_text("".join(f"{chrom[3:]}\t{p}\t{p}\n" for p in positions))
        for attempt in range(4):
            r = subprocess.run(["tabix", "-R", str(regions), url], capture_output=True, text=True, cwd=DATA)
            if r.returncode == 0:
                break
        else:
            raise RuntimeError(f"CADD tabix failed for {chrom}: {r.stderr[-300:]}")
        found = {}
        for line in r.stdout.splitlines():
            f = line.split("\t")
            key = ("chr" + f[0], int(f[1]), f[2], f[3])
            if key not in wanted:
                continue
            spl = [float(f[i]) for i in splice if f[i] != "NA"]
            prev = found.get(key)
            spl_max = max(spl) if spl else None
            if prev is None:
                found[key] = [f[col["PHRED"]], f[col["Consequence"]], f[col["GeneName"]], spl_max]
                continue
            # inclAnno has one row per annotation; keep the most severe consequence and the largest SpliceAI delta.
            if SEVERITY.get(f[col["Consequence"]], 99) < SEVERITY.get(prev[1], 99):
                prev[1], prev[2] = f[col["Consequence"]], f[col["GeneName"]]
            if spl_max is not None and (prev[3] is None or spl_max > prev[3]):
                prev[3] = spl_max
        print(f"  {chrom}: {len(found)}/{sum(1 for c, *_ in wanted if c == chrom)}", flush=True)
        return found

    rows = {}
    with ThreadPoolExecutor(6) as pool:
        for found in pool.map(one_chrom, sorted({c for c, *_ in wanted})):
            rows.update(found)
    print(f"  CADD {tag}: {len(rows)}/{len(wanted)} found")
    return rows


def revel(wanted: set) -> None:
    rows = {}
    with zipfile.ZipFile(REVEL_ZIP) as z, z.open("revel_with_transcript_ids") as raw:
        for line in io.TextIOWrapper(raw):
            f = line.split(",")
            if f[2] == "." or not f[2].isdigit():
                continue
            key = ("chr" + f[0], int(f[2]), f[3], f[4])
            if key in wanted:
                s = float(f[7])
                rows[key] = max(s, rows.get(key, s))
    write("revel", ["revel"], {k: [v] for k, v in rows.items()})


def alphamissense(wanted: set) -> None:
    rows = {}
    with gzip.open(AM, "rt") as h:
        for line in h:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]), f[2], f[3])
            if key in wanted:
                s = float(f[8])
                rows[key] = max(s, rows.get(key, s))
    write("am", ["alphamissense"], {k: [v] for k, v in rows.items()})


def gpn_msa(wanted: set) -> None:
    """tabix follows Hugging Face's redirect for the data but not for the index, so the index is fetched first."""
    idx = DATA / "gpn_msa.tsv.bgz.tbi"
    if not idx.exists():
        subprocess.run(["curl", "-sSfL", "-o", str(idx), GPN_MSA + ".tbi"], check=True)

    def one_chrom(chrom: str) -> dict:
        regions = DATA / f"kgp/gpn_msa_regions_{chrom}.tsv"
        regions.write_text("".join(f"{chrom[3:]}\t{p}\t{p}\n" for p in sorted({p for c, p, _, _ in wanted if c == chrom})))
        for attempt in range(6):
            r = subprocess.run(["tabix", "-R", str(regions), f"{GPN_MSA}##idx##{idx}"], capture_output=True, text=True)
            if r.returncode == 0:
                break
        else:
            raise RuntimeError(f"GPN-MSA tabix failed for {chrom}: {r.stderr[-300:]}")
        found = {}
        for line in r.stdout.splitlines():
            f = line.split("\t")
            key = ("chr" + f[0], int(f[1]), f[2], f[3])
            if key in wanted:
                found[key] = [-float(f[4])]
        print(f"  {chrom}: {len(found)}", flush=True)
        return found

    rows = {}
    with ThreadPoolExecutor(4) as pool:
        for found in pool.map(one_chrom, sorted({c for c, *_ in wanted})):
            rows.update(found)
    write("gpnmsa", ["gpn_msa"], rows)


def gpn_star(wanted: set) -> None:
    import duckdb
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs;")
    for model, name in (("v100", "gpnstar_v"), ("m447", "gpnstar_m")):
        rows = {}
        for chrom in sorted({c for c, *_ in wanted}):
            positions = sorted({p for c, p, _, _ in wanted if c == chrom})
            q = (f"SELECT pos, ref, alt, llr_calibrated FROM read_parquet('{GPN_STAR.format(model=model, chrom=chrom)}') "
                 f"WHERE pos IN ({','.join(map(str, positions))})")
            for pos, ref, alt, llr in con.sql(q).fetchall():
                key = (chrom, pos, ref, alt)
                if key in wanted and llr is not None:
                    rows[key] = [-llr]
        write(name, [name], rows)


def avi(wanted: set) -> None:
    snvs = {k for k in wanted if len(k[2]) == len(k[3]) == 1}
    regions = DATA / "kgp/avi_regions.tsv"
    regions.write_text("".join(f"{c}\t{p}\t{p}\n" for c, p in sorted({(c, p) for c, p, _, _ in snvs})))
    r = subprocess.run(["tabix", "-R", str(regions), str(AVI)], capture_output=True, text=True, check=True)
    rows = {}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        key = (f[0], int(f[1]), f[2], f[3])
        if key in snvs:
            rows[key] = [f[4], f[5]]
    print(f"  AVI: {len(rows)}/{len(snvs)} SNVs found")
    write("avi", ["avi", "avi_phred"], rows)


def spliceai_run(wanted: set) -> None:
    rows = {}
    with (DATA / "spliceai/indels.out.vcf").open() as h:
        for line in h:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]), f[3], f[4])
            if key not in wanted or not f[7].startswith("SpliceAI="):
                continue
            # One entry per gene: ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL; keep the largest delta.
            deltas = [float(x) for e in f[7][len("SpliceAI="):].split(",") for x in e.split("|")[2:6] if x != "."]
            if deltas:
                rows[key] = [max(deltas)]
    write("spliceai_run", ["spliceai_run"], rows)


def main() -> None:
    global OUT, REUSE
    p = argparse.ArgumentParser()
    p.add_argument("sources", nargs="*", default=["cadd", "revel", "am", "gpnmsa", "gpnstar"])
    p.add_argument("--study", default="snv", help="snv (data/) or a subdirectory of data/, e.g. indel")
    p.add_argument("--reuse", action="store_true", help="copy scores for variants in data/case_variants.tsv")
    a = p.parse_args()
    OUT = DATA if a.study == "snv" else DATA / a.study
    wanted = variants(OUT)
    print(f"{len(wanted)} case variants")
    if a.reuse and OUT != DATA:
        REUSE = wanted & variants(DATA)
        wanted -= REUSE
        print(f"{len(REUSE)} reused from the SNV study, {len(wanted)} to look up")
    which = a.sources
    for name in which:
        {"cadd": cadd, "revel": revel, "am": alphamissense, "gpnmsa": gpn_msa, "gpnstar": gpn_star,
         "spliceai_run": spliceai_run, "avi": avi}[name](wanted)


if __name__ == "__main__":
    main()
