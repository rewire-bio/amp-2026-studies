"""Mode of inheritance for every MANE gene in the panel regions, from PanelApp 285 green genes.

Implements PROTOCOL.md section 3. Each MANE gene symbol in region_genes.tsv is matched to a PanelApp 285 green gene
(confidence_level 3): first by Ensembl gene ID (PanelApp gene_data.ensembl_genes.GRch38.*.ensembl_id against the
MANE GTF gene_id, version stripped), then by symbol (PanelApp gene_symbol, hgnc_symbol, then alias). The PanelApp
mode_of_inheritance string is mapped to one category:

  mono      MONOALLELIC (any imprinting status) and BOTH monoallelic and biallelic
  biallelic BIALLELIC
  xl_mono   X-LINKED, "monoallelic mutations in females may cause disease"
  xl_bi     X-LINKED, "biallelic mutations in females"
  unknown   anything else, and genes with no green PanelApp match

Region rows that name two overlapping genes ("AKT3,SDCCAG8") are split, so each MANE symbol gets one row; callers
split GeneIndex names on commas too. Output data/gene_moi.tsv: gene (MANE symbol), ensg, moi, panelapp_moi, plus panelapp_gene and match (ensg, symbol,
alias or none) for checking.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python moi.py

New in the second run (no first-run counterpart).
"""
from __future__ import annotations

import csv
import gzip
import json
import re
from collections import Counter

import common

OUT = common.DATA / "gene_moi.tsv"


def category(moi: str) -> str:
    if moi.startswith("BIALLELIC"):
        return "biallelic"
    if moi.startswith("MONOALLELIC") or moi.startswith("BOTH"):
        return "mono"
    if moi.startswith("X-LINKED"):
        if "monoallelic mutations in females may cause disease" in moi:
            return "xl_mono"
        if "biallelic mutations in females" in moi:
            return "xl_bi"
    return "unknown"


def mane_gene_ids() -> dict[str, set[str]]:
    """MANE gene symbol to its Ensembl gene IDs (version stripped), from the gene lines of the MANE GTF."""
    out: dict[str, set[str]] = {}
    with gzip.open(common.MANE_GTF, "rt") as h:
        for line in h:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^".]+)', f[8]).group(1)
            name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
            out.setdefault(name, set()).add(gid)
    return out


def main() -> None:
    green = [g for g in json.loads(common.PANELAPP.read_text())["genes"]
             if str(g.get("confidence_level")) == "3" and g.get("entity_type") == "gene"]
    by_ensg: dict[str, list[dict]] = {}
    by_symbol: dict[str, list[dict]] = {}
    by_alias: dict[str, list[dict]] = {}
    for g in green:
        gd = g["gene_data"]
        for rel in (gd.get("ensembl_genes") or {}).get("GRch38", {}).values():
            if rel.get("ensembl_id"):
                by_ensg.setdefault(rel["ensembl_id"], []).append(g)
        for s in {gd.get("gene_symbol"), gd.get("hgnc_symbol"), g["entity_name"]} - {None, ""}:
            by_symbol.setdefault(s, []).append(g)
        for a in gd.get("alias") or []:
            by_alias.setdefault(a, []).append(g)

    # Where two MANE genes overlap, region_genes.tsv names both, comma-separated; each gets its own row.
    region_genes = sorted({g for r in csv.reader(common.REGION_GENES.open(), delimiter="\t") for g in r[3].split(",")})
    ensg_of = mane_gene_ids()
    rows, unmatched, multi, used = [], [], [], Counter()
    for gene in region_genes:
        ensgs = sorted(ensg_of.get(gene, set()))
        hits, how = [], "none"
        for e in ensgs:
            hits += by_ensg.get(e, [])
        if hits:
            how = "ensg"
        elif by_symbol.get(gene):
            hits, how = by_symbol[gene], "symbol"
        elif by_alias.get(gene):
            hits, how = by_alias[gene], "alias"
        names = sorted({h["entity_name"] for h in hits})
        if len(names) > 1:
            multi.append((gene, names))
        if hits:
            h = sorted(hits, key=lambda x: x["entity_name"])[0]
            used[h["entity_name"]] += 1
            rows.append([gene, ";".join(ensgs), category(h["mode_of_inheritance"]), h["mode_of_inheritance"],
                         h["entity_name"], how])
        else:
            unmatched.append(gene)
            rows.append([gene, ";".join(ensgs), "unknown", "", "", "none"])

    with OUT.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "ensg", "moi", "panelapp_moi", "panelapp_gene", "match"])
        w.writerows(rows)

    print(f"{len(region_genes)} MANE region genes, {len(green)} PanelApp green genes")
    print("matched by:", dict(Counter(r[5] for r in rows)))
    print("moi:", dict(Counter(r[2] for r in rows)))
    print("no MANE gene ID:", [g for g in region_genes if g not in ensg_of])
    print("unmatched:", unmatched)
    print("several PanelApp genes for one MANE gene:", multi)
    print("PanelApp genes used by two MANE genes:", [g for g, n in used.items() if n > 1])
    print("green genes not matched to any region gene:", sorted({g["entity_name"] for g in green} - set(used)))
    for r in rows:
        if r[5] in ("symbol", "alias") or r[0] != r[4] and r[4]:
            print("  check:", r[0], r[1], "->", r[4], f"({r[5]})")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
