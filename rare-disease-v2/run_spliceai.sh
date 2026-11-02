#!/usr/bin/env bash
# Run SpliceAI 1.3.1 on the indels of the indel study, as labs do for indels missing from the precomputed files
# (PROTOCOL.md section 6). CPU only, about an hour per 2,000 new indels. Run from rare-disease-v2 after
# build_cases.py, then lookup_scores.py --study indel spliceai_run.
# Derived from rare-disease/run_spliceai.sh at amp-2026-v4.
#
#   bash run_spliceai.sh
#
# Input: the indels (REF or ALT longer than 1 bp) in data/indel/case_variants.tsv. Output: data/spliceai/indels.out.vcf.
# Needs an uncompressed, indexed hg38 FASTA (HG38, default data/spliceai/hg38.fa; see PLAN.md phase 3).
#
# Reuse: the first run scored its indels with the same SpliceAI version, the same command and the same hg38
# (../rare-disease/data/spliceai/indels.out.vcf, which has one row per input indel, INFO "." where SpliceAI gives
# no score). Rows for indels already in that file, or in this run's own earlier output ($OUT, made the same way),
# are copied, and SpliceAI runs only on the new ones. REUSE=0 runs SpliceAI on every indel.
#
# Overrides for tests: VARIANTS (any table with chrom, pos, ref, alt), OUT (output VCF), HG38, REUSE (1 or 0).
set -euo pipefail
cd "$(dirname "$0")"
VARIANTS=${VARIANTS:-data/indel/case_variants.tsv}
OUT=${OUT:-data/spliceai/indels.out.vcf}
HG38=${HG38:-data/spliceai/hg38.fa}
REUSE=${REUSE:-1}
V1_OUT=../rare-disease/data/spliceai/indels.out.vcf
[ -s "$HG38.fai" ] || { echo "no FASTA index $HG38.fai: decompress and samtools faidx hg38 first" >&2; exit 1; }
WORK="$OUT.work"
mkdir -p "$WORK"

# Split the indels into those already scored (rows copied) and new ones (written to new.vcf for SpliceAI).
python3 - "$VARIANTS" "$V1_OUT" "$REUSE" "$WORK" "$OUT" <<'PY'
import csv, os, sys
variants, v1_out, reuse, work, own_out = sys.argv[1], sys.argv[2], sys.argv[3] == "1", sys.argv[4], sys.argv[5]
keys = set()
with open(variants) as h:
    for r in csv.DictReader(h, delimiter="\t"):
        if not (len(r["ref"]) == 1 and len(r["alt"]) == 1):
            keys.add((r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
reused, source = {}, {}
if reuse:
    # This run's earlier output first, then the first run's, so a row is taken from this run where both have it.
    for name, path in (("this run", own_out), ("first run", v1_out)):
        if not os.path.exists(path):
            continue
        with open(path) as h:
            for line in h:
                if not line.startswith("#"):
                    f = line.split("\t")
                    k = (f[0], int(f[1]), f[3], f[4])
                    if k in keys and k not in reused:
                        reused[k], source[k] = line, name
order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
srt = lambda k: (order[k[0]], k[1], k[2], k[3])
hdr = "##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={c}>\n" for c in order) + \
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
with open(f"{work}/new.vcf", "w") as out:
    out.write(hdr)
    for c, p, r, a in sorted(keys - set(reused), key=srt):
        out.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n")
with open(f"{work}/reused.body", "w") as out:
    out.writelines(reused[k] for k in sorted(reused, key=srt))
n_own = sum(1 for v in source.values() if v == "this run")
print(f"{len(keys)} indels: {n_own} copied from this run's earlier output, {len(reused) - n_own} from the first run, "
      f"{len(keys) - len(reused)} to score")
PY

if grep -qv '^#' "$WORK/new.vcf"; then
  # SpliceAI 1.3.1 needs NumPy 1 (np.fromstring) and Keras 2 (TensorFlow 2.15).
  uv run -q --python 3.11 --with spliceai==1.3.1 --with "tensorflow==2.15.*" --with "numpy<2" \
    spliceai -I "$WORK/new.vcf" -O "$WORK/new.out.vcf" -R "$HG38" -A grch38
  HEADER="$WORK/new.out.vcf"
else
  HEADER="$V1_OUT"
  : > "$WORK/new.out.vcf"
fi

# One header (SpliceAI's), then every row in contig order.
python3 - "$HEADER" "$WORK/new.out.vcf" "$WORK/reused.body" "$OUT" <<'PY'
import sys
header_src, new_out, reused, out = sys.argv[1:]
order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
header = [l for l in open(header_src) if l.startswith("#")]
body = [l for l in open(new_out) if not l.startswith("#")] + open(reused).readlines()
body.sort(key=lambda l: (order[l.split("\t")[0]], int(l.split("\t")[1]), l.split("\t")[3], l.split("\t")[4]))
with open(out, "w") as h:
    h.writelines(header + body)
print(f"wrote {out}: {len(body)} indels")
PY
rm -r "$WORK"
