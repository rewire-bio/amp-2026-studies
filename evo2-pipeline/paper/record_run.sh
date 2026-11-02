#!/usr/bin/env bash
# Run the pipeline's replay profile in a scratch copy and record what the paper reports in run/.
#
#   bash evo2-pipeline/paper/record_run.sh [path to evo2-cloud-operations]   # default ../../evo2-cloud-operations
#
# Runs, each in a fresh work directory unless noted:
#   first, second    replay with the H100 (GCP) scores, twice, to compare outputs
#   resume           the second run again with -resume, to count cached tasks
#   p3               replay with the same variants' scores from the Runpod H100
#   l40s             replay with the same variants' scores from the Runpod L40S
#   mismatch         the test case plus one record whose REF does not match the genome
#   missing          the test case plus one SNV with no recorded Evo 2 score (the run must stop)
# After the first run, the windows MAKE_WINDOWS built are compared with the windows evo2-cloud-operations builds for
# the same variants (src/evo2ops/genome.py snv_windows, from that study's raw chromosome 17 FASTA), by SHA-256.
# then summarise_run.py writes run/summary.json, run/ranks.tsv and run/trace_first.txt. CADD and SpliceAI values are
# not written to run/ (their licences do not allow redistribution); their ranks are.
#
# Needs Nextflow 24.10 or later with Java 17+ (JAVA_HOME defaults to Homebrew's openjdk@21), bcftools, tabix, bgzip,
# samtools and python3. Nothing in the study folder is changed except run/.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIPE="$(dirname "$HERE")"
OPS="$(cd "${1:-$PIPE/../../evo2-cloud-operations}" && pwd)"
export JAVA_HOME="${JAVA_HOME:-/opt/homebrew/opt/openjdk@21}"
export NXF_ANSI_LOG=false
REC="$OPS/evidence/runs/recorded-b75c67aa/platforms"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

mkdir -p "$W/p"
cp -R "$PIPE/main.nf" "$PIPE/nextflow.config" "$PIPE/bin" "$PIPE/scripts" "$PIPE/test" "$W/p/"
(cd "$W/p" && python3 scripts/prepare_resources.py --study "$OPS")
# prepare_resources.py rewrites test/ from its fixed seed: it must match the committed case
if diff -r "$W/p/test" "$PIPE/test" >/dev/null; then echo identical > "$W/test_check.txt"; else echo differs > "$W/test_check.txt"; fi

# Two extra cases for the checks: one REF mismatch, one SNV without a recorded score. 17:41250000 is in a BRCA1
# intron and not among the assay's variants.
BASE="$(samtools faidx "$W/p/resources/GRCh37_chr17.fa.gz" 17:41250000-41250000 | tail -1 | tr a-z A-Z)"
python3 - "$W/p/test/case_brca1.vcf" "$BASE" "$W" <<'PY'
import sys
case, base, out = sys.argv[1], sys.argv[2], sys.argv[3]
lines = open(case).read().splitlines()
head = [l for l in lines if l.startswith("#")]
body = [l for l in lines if not l.startswith("#")]
other = next(b for b in "ACGT" if b != base)
alt = next(b for b in "ACGT" if b not in (base, other))
def write(name, extra):
    rows = sorted(body + [extra], key=lambda l: int(l.split("\t")[1]))
    open(f"{out}/{name}.vcf", "w").write("\n".join(head + rows) + "\n")
write("mismatch", f"17\t41250000\t.\t{other}\t{alt}\t.\tPASS\t.")
write("missing", f"17\t41250000\t.\t{base}\t{other}\t.\tPASS\t.")
PY

run() {  # name, then extra arguments for nextflow run
  local name="$1"; shift
  local start end status=0
  start=$(python3 -c 'import time; print(time.time())')
  (cd "$W/p" && nextflow run main.nf -profile replay -with-trace "$W/trace_$name.txt" "$@" > "$W/log_$name.txt" 2>&1) || status=$?
  end=$(python3 -c 'import time; print(time.time())')
  python3 -c "print(round($end - $start, 2))" > "$W/wall_$name.txt"
  echo "$status" > "$W/exit_$name.txt"
  mkdir -p "$W/out_$name"
  cp "$W/p/results/"* "$W/out_$name/" 2>/dev/null || true
  # keep the step outputs the checks read
  for f in normalised.vcf.gz windows.tsv evo2_scores.tsv; do
    find "$W/p/work" -name "$f" -exec cp {} "$W/out_$name/" \; 2>/dev/null || true
  done
  find "$W/p/work" -path '*/.command.err' -exec cat {} \; > "$W/err_$name.txt" 2>/dev/null || true
  echo "$name: exit $status"
}
fresh() { rm -rf "$W/p/work" "$W/p/results" "$W/p/.nextflow" "$W/p"/.nextflow.log*; }

fresh; run first
# The windows the pipeline built must be the sequences the recorded GPU run scored: rebuild them with the source
# study's own code and compare.
python3 - "$OPS" "$W/out_first/windows.tsv" "$W/windows_check.json" <<'WIN'
import csv, hashlib, json, sys
from pathlib import Path
ops, pipeline_windows, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
sys.path.insert(0, str(ops / "src"))
from evo2ops.genome import read_first_fasta_record, snv_windows
chromosome = read_first_fasta_record(ops / "data/raw/GRCh37.p13_chr17.fna.gz")
with open(pipeline_windows) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
def digest(pairs):
    sequences = sorted({s for pair in pairs for s in pair})
    return hashlib.sha256("\n".join(sequences).encode()).hexdigest(), len(sequences)
mine = [(r["ref_seq"], r["var_seq"]) for r in rows]
theirs = [snv_windows(chromosome, int(r["pos"]), r["ref"], r["alt"], 8192) for r in rows]
(d1, n1), (d2, n2) = digest(mine), digest(theirs)
same = sum(a == b for a, b in zip(mine, theirs))
json.dump({"variants": len(rows), "windows_pipeline": n1, "windows_source": n2, "variants_identical": same,
           "sha256_pipeline": d1, "sha256_source": d2, "identical": d1 == d2 and same == len(rows)},
          open(out, "w"), indent=1)
WIN
fresh; run second
run resume -resume
fresh; run p3 --evo2_recorded "$REC/P3/attempt-1/w1_brca1_scores.tsv"
fresh; run l40s --evo2_recorded "$REC/P5/attempt-2/w1_brca1_scores.tsv"
fresh; run mismatch --case_vcf "$W/mismatch.vcf"
fresh; run missing --case_vcf "$W/missing.vcf"

{
  echo "nextflow	$(cd "$W/p" && nextflow -version 2>&1 | sed -nE 's/.*version ([0-9.]+) build ([0-9]+).*/\1 build \2/p')"
  echo "java	$("$JAVA_HOME/bin/java" -version 2>&1 | head -1)"
  echo "bcftools	$(bcftools --version | head -1)"
  echo "htslib	$(bcftools --version | sed -n 2p)"
  echo "samtools	$(samtools --version | head -1)"
  echo "python	$(python3 --version)"
  echo "os	$(uname -sr) $(uname -m)"
  echo "evo2-cloud-operations	$(git -C "$OPS" rev-parse HEAD)"
  echo "amp-2026-studies	$(git -C "$PIPE" rev-parse HEAD)"
} > "$W/versions.tsv"

python3 "$HERE/summarise_run.py" --work "$W" --ops "$OPS" --out "$HERE/run"
