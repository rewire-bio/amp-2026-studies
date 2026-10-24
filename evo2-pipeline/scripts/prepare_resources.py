"""Build the pipeline's local resources and test case from rewire-bio/evo2-cloud-operations.

  python3 scripts/prepare_resources.py --study ../../evo2-cloud-operations

Writes resources/ (not committed: CADD and the assay labels are licensed for non-commercial use only):
  GRCh37_chr17.fa.gz (+ .fai, .gzi)   chr17 renamed to "17", bgzipped and indexed
  cadd_v1.7.tsv.gz (+ .tbi)           CADD v1.7 raw scores for the BRCA1 region, same layout as CADD's file
  spliceai_max.tsv.gz (+ .tbi)        largest of the four SpliceAI masked deltas, BRCA1 transcripts
  evo2_7b_h100_brca1.tsv              recorded Evo 2 7B scores from the study's H100 run
and test/: one simulated case, a loss-of-function variant plus 99 other assay variants at other positions,
with the diagnostic variant named in test/case_brca1.truth.txt. Needs bgzip, tabix and samtools.
"""
import argparse
import csv
import gzip
import json
import math
import random
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
RES, TEST = HERE / "resources", HERE / "test"
SPLICEAI = ["SpliceAI-acc-gain", "SpliceAI-acc-loss", "SpliceAI-don-gain", "SpliceAI-don-loss"]
SEED, CANDIDATES = 20261016, 99


def bgzip_tabix(rows: list[list[str]], header: list[str], out: Path) -> None:
    rows.sort(key=lambda r: (r[0], int(r[1])))
    plain = out.with_suffix("")
    with plain.open("w") as handle:
        handle.write("#" + "\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", str(out)], check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", type=Path, required=True)
    study = parser.parse_args().study.resolve()
    raw = study / "data/raw/current-predictors"
    RES.mkdir(exist_ok=True)
    TEST.mkdir(exist_ok=True)

    fasta = RES / "GRCh37_chr17.fa"
    with gzip.open(study / "data/raw/GRCh37.p13_chr17.fna.gz", "rt") as src, fasta.open("w") as dst:
        for line in src:
            dst.write(">17\n" if line.startswith(">") else line)
    subprocess.run(["bgzip", "-f", str(fasta)], check=True)
    subprocess.run(["samtools", "faidx", str(fasta) + ".gz"], check=True)

    with (raw / "cadd_v1.7.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cadd = [[r["Chrom"], r["Pos"], r["Ref"], r["Alt"], r["RawScore"]] for r in reader]
    bgzip_tabix(cadd, ["Chrom", "Pos", "Ref", "Alt", "RawScore"], RES / "cadd_v1.7.tsv.gz")

    best: dict[tuple, float] = {}
    with (raw / "spliceai_via_cadd_v1.7.tsv").open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            values = [float(r[c]) for c in SPLICEAI if r[c] not in ("NA", "")]
            if r["GeneName"] == "BRCA1" and values:
                key = (r["Chrom"], r["Pos"], r["Ref"], r["Alt"])
                best[key] = max(best.get(key, -math.inf), max(values))
    bgzip_tabix([[*k, repr(v)] for k, v in best.items()], ["Chrom", "Pos", "Ref", "Alt", "SpliceAI_max"],
                RES / "spliceai_max.tsv.gz")

    run = study / "evidence/runs/recorded-b75c67aa"
    attempt = json.loads((run / "operations.json").read_text())["platforms"]["P2"]["attempt"]
    shutil.copy(run / "platforms/P2" / attempt / "w1_brca1_scores.tsv", RES / "evo2_7b_h100_brca1.tsv")

    with (study / "data/derived/brca1_variants.tsv").open() as handle:
        variants = list(csv.DictReader(handle, delimiter="\t"))
    rng = random.Random(SEED)
    diagnostic = rng.choice([v for v in variants if v["lof"] == "1"])
    others = [v for v in variants if v["lof"] == "0" and v["pos"] != diagnostic["pos"]]
    case = sorted([diagnostic, *rng.sample(others, CANDIDATES)], key=lambda v: int(v["pos"]))
    with (TEST / "case_brca1.vcf").open("w") as out:
        out.write("##fileformat=VCFv4.2\n##reference=GRCh37\n##contig=<ID=17,length=81195210>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in case:
            out.write(f"17\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}\t.\tPASS\t.\n")
    (TEST / "case_brca1.truth.txt").write_text(f"17:{diagnostic['pos']}:{diagnostic['ref']}>{diagnostic['alt']}\t{diagnostic['consequence']}\n")
    print(f"resources in {RES}; test case of {len(case)} variants, diagnostic {diagnostic['id']} ({diagnostic['consequence']})")


if __name__ == "__main__":
    main()
