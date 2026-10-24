#!/usr/bin/env bash
# Run SpliceAI 1.3.1 on every indel in the indel and cancer studies, as labs do for indels missing from the
# precomputed files. CPU only, about an hour for ~2,200 indels. Run from rare-disease after build_cases.py.
# Needs an uncompressed, indexed hg38 FASTA (HG38, default data/spliceai/hg38.fa).
set -euo pipefail
cd "$(dirname "$0")/data"
mkdir -p spliceai
HG38=${HG38:-spliceai/hg38.fa}
python3 - <<'PY'
import csv
keys = set()
for st in ("indel", "cancer"):
    for r in csv.DictReader(open(f"{st}/case_variants.tsv"), delimiter="\t"):
        if not (len(r["ref"]) == 1 and len(r["alt"]) == 1):
            keys.add((r["chrom"], int(r["pos"]), r["ref"], r["alt"]))
order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
hdr = "##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={c}>\n" for c in order) + \
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
with open("spliceai/indels.vcf", "w") as out:
    out.write(hdr)
    for c, p, r, a in sorted(keys, key=lambda k: (order[k[0]], k[1], k[2], k[3])):
        out.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n")
print(len(keys), "indels")
PY
# SpliceAI 1.3.1 needs NumPy 1 (np.fromstring) and Keras 2 (TensorFlow 2.15).
uv run -q --python 3.11 --with spliceai==1.3.1 --with "tensorflow==2.15.*" --with "numpy<2" \
  spliceai -I spliceai/indels.vcf -O spliceai/indels.out.vcf -R "$HG38" -A grch38
