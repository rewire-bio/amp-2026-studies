"""Target regions: MANE Select exons (plus 20 bp each side) of a PanelApp panel's green genes.

  python3 rare-disease/regions.py                  # panel 285 Intellectual disability, writes data/regions.tsv
  python3 rare-disease/regions.py 245 cancer       # panel 245, writes data/cancer/regions.tsv

regions.tsv is merged and 1-based inclusive, for bcftools -R. The panel JSON is read from the output directory as
panelapp_<id>.json.
"""
import gzip
import json
import re
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
PANEL = sys.argv[1] if len(sys.argv) > 1 else "285"
D = DATA / sys.argv[2] if len(sys.argv) > 2 else DATA
PAD = 20
panel = json.loads((D / f"panelapp_{PANEL}.json").read_text())
green, green_ids = set(), set()
for g in panel["genes"]:
    if g["confidence_level"] != "3":
        continue
    green.add(g["entity_name"])
    for build in ("GRch38", "GRCh38"):
        for rel in (g["gene_data"].get("ensembl_genes") or {}).get(build, {}).values():
            green_ids.add(rel["ensembl_id"])
by_chrom = {}
genes_found = set()
with gzip.open(DATA / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz", "rt") as h:
    for line in h:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or 'tag "MANE_Select"' not in f[8]:
            continue
        gene = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
        gene_id = re.search(r'gene_id "(ENSG\d+)', f[8]).group(1)
        if gene not in green and gene_id not in green_ids:
            continue
        genes_found.add(gene)
        by_chrom.setdefault(f[0], []).append((max(1, int(f[3]) - PAD), int(f[4]) + PAD, gene))
order = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
n = 0
with (D / "regions.tsv").open("w") as out, (D / "region_genes.tsv").open("w") as genes_out:
    for chrom in order:
        merged = []
        for s, e, g in sorted(by_chrom.get(chrom, [])):
            if merged and s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
                merged[-1][2].add(g)
            else:
                merged.append([s, e, {g}])
        for s, e, gs in merged:
            out.write(f"{chrom}\t{s}\t{e}\n")
            genes_out.write(f"{chrom}\t{s}\t{e}\t{','.join(sorted(gs))}\n")
            n += 1
print(f"{len(green)} green genes, {len(genes_found)} matched to a MANE Select transcript, {n} merged regions")
