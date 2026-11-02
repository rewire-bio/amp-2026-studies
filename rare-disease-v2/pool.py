"""ClinVar diagnostic candidate pools: class, frequency, gene and inheritance mode for every candidate.

Implements PROTOCOL.md sections 3, 4 and 5 for the pools (not the case build). Inputs are the first run's ClinVar
P/LP 2+ star tables, ../rare-disease/data/clinvar_plp2.tsv (SNVs) and clinvar_indel_plp2.tsv (indels up to 50 bp).
The indel pool leaves out variants with REF and ALT both longer than 1 bp (multi-base substitutions and complex
changes, section 5); how many, by class, is logged and written to results/pool_counts.md.

- SNV class (section 4): CADD v1.7's most severe consequence over MANE Select rows, else over all rows
  (common.cadd_lookup). Cached in data/pool/scores_cadd.tsv (gitignored; CADD may not be redistributed).
- Indel class: annotate_indels.py on the indel pool, written to data/pool/consequence_indels.tsv.
- Full gnomAD v3.1.2 genome counts (common.myvariant_counts) for both pools, cached in data/pool/gnomad_full.tsv in
  the layout of data/gnomad_full.tsv. A variant not found has grpmax 0. Eligible: grpmax below common.MAIN_AF.
- Gene: the MANE gene of the region (common.GeneIndex; names of overlapping genes are split on commas). Where
  several genes match, the CADD MANE gene is preferred for SNVs, then the ClinVar gene, then the first by name.
  Inheritance mode from data/gene_moi.tsv (moi.py); "none" where the variant is in no region.

Class names (an interface for build_cases.py and case_benchmark.py):
  SNV:   nonsense, missense, canonical_splice, splice_region, synonymous, utr, intronic, other
  indel: frameshift, inframe, canonical_splice (annotate_indels "splice"), utr, intronic, other
Diagnostic groups:
  SNV:   missense, canonical_splice, nonsense, noncoding (splice_region, synonymous, utr, intronic)
  indel: frameshift, inframe, canonical_splice, noncoding (utr, intronic)
  Anything else has group "none".

Outputs: data/pool/snv_pool.tsv and data/pool/indel_pool.tsv (the ClinVar columns, with ClinVar's gene renamed
clinvar_gene, plus class, group, cadd_consequence and cadd_mane for SNVs, grpmax, found, gene, moi, eligible), and
results/pool_counts.md (eligible variants and positions per group against the quotas of section 5).

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python pool.py

New in the second run. The SNV classes replace rare-disease/build_cases.py's group() on ClinVar's MC field
(amp-2026-v4); the frequency lookup follows rare-disease/gnomad_af.py at amp-2026-v4 through common.py.
"""
from __future__ import annotations

import argparse
import csv
import random
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import annotate_indels
import common

POOL = common.DATA / "pool"
SNV_IN = common.V1_DATA / "clinvar_plp2.tsv"
INDEL_IN = common.V1_DATA / "clinvar_indel_plp2.tsv"
CADD_CACHE = POOL / "scores_cadd.tsv"
GNOMAD_CACHE = POOL / "gnomad_full.tsv"
INDEL_CLASSES = POOL / "consequence_indels.tsv"
FIRST_RUN_CASES = [common.V1_DATA / "cases.tsv", common.V1_DATA / "replication" / "cases.tsv"]
COUNTS_MD = common.RESULTS / "pool_counts.md"

SNV_CLASS = {"STOP_GAINED": "nonsense", "NON_SYNONYMOUS": "missense", "CANONICAL_SPLICE": "canonical_splice",
             "SPLICE_SITE": "splice_region", "SYNONYMOUS": "synonymous", "5PRIME_UTR": "utr", "3PRIME_UTR": "utr",
             "INTRONIC": "intronic"}
INDEL_CLASS = {"splice": "canonical_splice", "frameshift": "frameshift", "inframe": "inframe", "utr": "utr",
               "intronic": "intronic", "other": "other"}
SNV_GROUP = {"missense": "missense", "canonical_splice": "canonical_splice", "nonsense": "nonsense",
             "splice_region": "noncoding", "synonymous": "noncoding", "utr": "noncoding", "intronic": "noncoding"}
INDEL_GROUP = {"frameshift": "frameshift", "inframe": "inframe", "canonical_splice": "canonical_splice",
               "utr": "noncoding", "intronic": "noncoding"}
# Section 5 quotas. SNV: the 200 and the 500 draw from the same pool.
Q200 = {"missense": 60, "canonical_splice": 40, "nonsense": 30, "noncoding": 70}
Q500 = {"missense": 150, "canonical_splice": 100, "nonsense": 75, "noncoding": 175}
Q120 = {"frameshift": 40, "inframe": 30, "canonical_splice": 30, "noncoding": 20}
GENE_CAP = 10
CADD_COLS = ["cadd_phred", "cadd_consequence", "cadd_consequence_any", "gene", "mane", "spliceai_max"]


def log(*a) -> None:
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def key_of(r: dict) -> tuple:
    return (r["chrom"], int(r["pos"]), r["ref"], r["alt"])


def write_tsv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


# ---------------------------------------------------------------------------------------------------------------
# Cached lookups

def cadd_chrom(chrom: str, keys: set[tuple], workdir: Path) -> dict[tuple, dict]:
    """common.cadd_lookup for one chromosome, with a pause between attempts and the result kept in workdir.

    The CADD server drops some connections when several remote tabix reads run at once, and cadd_lookup retries
    without waiting, so failures are retried here after a pause. A rerun reuses finished chromosomes.
    """
    import json
    part = workdir / f"pool_{chrom}.json"
    if part.exists():
        return {tuple(json.loads(k)): v for k, v in json.loads(part.read_text()).items()}
    for attempt in range(8):
        try:
            rows = common.cadd_lookup(keys, workdir, threads=1)
            break
        except RuntimeError as e:
            log(f"CADD {chrom}: attempt {attempt + 1} failed ({str(e)[:80]}), retrying")
            time.sleep(30 * (attempt + 1))
    else:
        raise SystemExit(f"CADD lookup failed for {chrom}")
    part.write_text(json.dumps({json.dumps(list(k)): v for k, v in rows.items()}))
    log(f"CADD {chrom}: {len(rows)} of {len(keys)} found")
    return rows


def cadd(snv_keys: list[tuple], workdir: Path) -> dict[tuple, dict]:
    if not CADD_CACHE.exists():
        from concurrent.futures import ThreadPoolExecutor
        log(f"CADD: looking up {len(snv_keys)} SNVs by remote tabix (workdir {workdir})")
        workdir.mkdir(parents=True, exist_ok=True)
        by_chrom = defaultdict(set)
        for k in snv_keys:
            by_chrom[k[0]].add(k)
        rows = {}
        with ThreadPoolExecutor(3) as ex:
            for part in ex.map(lambda c: cadd_chrom(c, by_chrom[c], workdir), sorted(by_chrom)):
                rows.update(part)
        write_tsv(CADD_CACHE, common.KEY + CADD_COLS,
                  ([*k, *(fmt(rows[k][c]) for c in CADD_COLS)] for k in snv_keys if k in rows))
        log(f"CADD: {len(rows)} of {len(snv_keys)} found, written to {CADD_CACHE}")
    return {key_of(r): r for r in read_tsv(CADD_CACHE)}


def gnomad(keys: list[tuple]) -> dict[tuple, dict]:
    cols = ["found", *common.COUNT_FIELDS, "af", "grpmax"]
    if not GNOMAD_CACHE.exists():
        log(f"gnomAD: querying MyVariant.info for {len(keys)} variants")
        counts = common.myvariant_counts(keys)
        out = []
        for k in keys:
            d = counts.get(k)
            if d is None:
                out.append([*k, "False", *(["0"] * len(common.COUNT_FIELDS)), "0", "0"])
            else:
                out.append([*k, "True", *(str(d[c]) for c in common.COUNT_FIELDS), fmt(d["af"]), fmt(d["grpmax"])])
        write_tsv(GNOMAD_CACHE, common.KEY + cols, out)
        log(f"gnomAD: {len(counts)} of {len(keys)} found, written to {GNOMAD_CACHE}")
    return {key_of(r): r for r in read_tsv(GNOMAD_CACHE)}


# ---------------------------------------------------------------------------------------------------------------
# Pools

# ClinVar genes whose disease is outside the panel but whose variants fall inside a panel gene's region (MSH3 inside
# DHFR, MUTYH inside TOE1). Their variants are not eligible: the case would be a diagnosis the panel does not cover.
OFF_PANEL_GENES = {"MSH3", "MUTYH"}


def off_panel(clinvar_gene: str) -> bool:
    genes = {g.split(":")[0] for g in clinvar_gene.split("|") if g}
    return bool(genes) and genes <= OFF_PANEL_GENES


def gene_for(r: dict, index: common.GeneIndex, prefer: list[str]) -> tuple[str, int]:
    """MANE gene for a pool row, and the number of candidate genes."""
    end = int(r["pos"]) + max(len(r["ref"]) - 1, 1 if len(r["alt"]) > 1 else 0)
    cands = sorted({g for name in index.genes(r["chrom"], int(r["pos"]), end) for g in name.split(",")})
    for p in prefer:
        if p in cands:
            return p, len(cands)
    return (cands[0] if cands else ""), len(cands)


def build(kind: str, rows: list[dict], classes: dict[tuple, str], cadd_rows: dict[tuple, dict],
          af: dict[tuple, dict], index: common.GeneIndex, moi: dict[str, str]) -> list[dict]:
    group_map = SNV_GROUP if kind == "snv" else INDEL_GROUP
    out, ngenes = [], Counter()
    for r in rows:
        k = key_of(r)
        c = cadd_rows.get(k, {})
        prefer = [c.get("gene", ""), r["gene"]] if kind == "snv" else [r["gene"]]
        gene, n = gene_for(r, index, [p for p in prefer if p])
        ngenes[min(n, 2)] += 1
        g = af[k]
        row = {**{("clinvar_gene" if f == "gene" else f): v for f, v in r.items()},
               "class": classes[k], "group": group_map.get(classes[k], "none")}
        if kind == "snv":
            row["cadd_consequence"] = c.get("cadd_consequence", "")
            row["cadd_mane"] = c.get("mane", "")
        row.update(grpmax=g["grpmax"], found=g["found"], gene=gene, moi=moi.get(gene, "none") if gene else "none",
                   eligible=int(float(g["grpmax"]) < common.MAIN_AF and not off_panel(r["gene"])))
        row["found"] = "1" if g["found"] in ("1", "True") else "0"
        if kind == "snv":
            row["cadd_mane"] = "1" if row["cadd_mane"] in ("1", "True") else ("0" if row["cadd_mane"] else "")
        out.append(row)
    log(f"{kind}: genes per variant: none {ngenes[0]}, one {ngenes[1]}, several {ngenes[2]}")
    return out


# ---------------------------------------------------------------------------------------------------------------
# Counts against the quotas

def first_run_positions() -> set[tuple]:
    pos = set()
    for f in FIRST_RUN_CASES:
        for r in read_tsv(f):
            if r["diagnostic"] == "1":
                pos.add((r["chrom"], int(r["pos"])))
    return pos


def positions_by_gene(rows: list[dict]) -> dict[tuple, str]:
    """Distinct positions of eligible rows and the gene of each (the first row's gene)."""
    out = {}
    for r in rows:
        out.setdefault((r["chrom"], int(r["pos"])), r["gene"])
    return out


def capped(per_gene: Counter) -> int:
    return sum(min(n, GENE_CAP) for n in per_gene.values())


def draw(pos_gene: dict[tuple, str], n: int, rng: random.Random) -> list[tuple]:
    """n positions at random with at most GENE_CAP per gene (as build_cases.py would draw)."""
    order = sorted(pos_gene)
    rng.shuffle(order)
    taken, per_gene = [], Counter()
    for p in order:
        if len(taken) == n:
            break
        if per_gene[pos_gene[p]] < GENE_CAP:
            taken.append(p)
            per_gene[pos_gene[p]] += 1
    return taken


def snv_section(pool: list[dict], first: set[tuple]) -> tuple[list[str], dict]:
    lines, verdict = [], {}
    elig = [r for r in pool if r["eligible"] == 1]
    lines += ["## SNV pool", "",
              f"{len(pool)} ClinVar SNVs; {len(elig)} eligible (full gnomAD grpmax below {common.MAIN_AF:.1%}). "
              f"First-run diagnostic positions: {len(first)}.", "",
              "The 200 may use any eligible position. The 500 may not use first-run diagnostic positions or this "
              "run's 200. Cap: at most 10 diagnostics per gene per group in each study. Capacity under the cap is "
              "the sum over genes of min(positions, 10).", "",
              "| Group | Eligible variants | Positions | Not first-run | Genes | Genes over 10 | Cap capacity, 200 "
              "| Cap capacity, 500 (not first-run) | Quota 200 | Quota 500 | Together |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    sims = {}
    for grp in Q200:
        rows = [r for r in elig if r["group"] == grp]
        pg = positions_by_gene(rows)
        new = {p: g for p, g in pg.items() if p not in first}
        genes, genes_new = Counter(pg.values()), Counter(new.values())
        over = sum(1 for n in genes.values() if n > GENE_CAP)
        lines.append(f"| {grp} | {len(rows)} | {len(pg)} | {len(new)} | {len(genes)} | {over} | {capped(genes)} "
                     f"| {capped(genes_new)} | {Q200[grp]} | {Q500[grp]} | {Q200[grp] + Q500[grp]} |")
        # Simulate the 200 then the 500, both capped, 200 times.
        short = []
        for seed in range(200):
            rng = random.Random(seed)
            used = set(draw(pg, Q200[grp], rng))
            left = {p: g for p, g in new.items() if p not in used}
            short.append(min(len(left), capped(Counter(left.values()))) - Q500[grp])
        sims[grp] = (min(short), sorted(short)[len(short) // 2], max(short))
        verdict[grp] = {"positions": len(pg), "new": len(new), "cap200": capped(genes),
                        "cap500": capped(genes_new), "over": over, "genes": len(genes),
                        "top": genes.most_common(3), "top_new": genes_new.most_common(3)}
    lines += ["", "Simulated draws (200 seeds): the 200 drawn at random under the cap, then the positions left for "
              "the 500 (not first-run, not the 200) under the cap, minus the 500's quota. Negative means a shortfall.",
              "", "| Group | Spare for the 500: min | median | max | Largest genes (positions) |", "|---|---|---|---|---|"]
    for grp, (lo, med, hi) in sims.items():
        top = ", ".join(f"{g} {n}" for g, n in verdict[grp]["top"])
        lines.append(f"| {grp} | {lo} | {med} | {hi} | {top} |")
        verdict[grp]["sim"] = (lo, med, hi)
    return lines, verdict


def indel_section(pool: list[dict], multi: Counter) -> tuple[list[str], dict]:
    lines, verdict = [], {}
    elig = [r for r in pool if r["eligible"] == 1]
    lines += ["## Indel pool", "",
              f"{len(pool)} ClinVar indels up to 50 bp, not counting {sum(multi.values())} variants with REF and ALT "
              f"both longer than 1 bp, which are left out (by class: "
              f"{', '.join(f'{c} {n}' for c, n in sorted(multi.items(), key=lambda x: INDEL_CLASS_ORDER.index(x[0])))}"
              f"); {len(elig)} eligible.", "",
              "| Group | Eligible variants | Positions | Genes | Genes over 10 | Cap capacity | Quota 120 | Spare |",
              "|---|---|---|---|---|---|---|---|"]
    for grp, q in Q120.items():
        rows = [r for r in elig if r["group"] == grp]
        pg = positions_by_gene(rows)
        genes = Counter(pg.values())
        over = sum(1 for n in genes.values() if n > GENE_CAP)
        cap = capped(genes)
        lines.append(f"| {grp} | {len(rows)} | {len(pg)} | {len(genes)} | {over} | {cap} | {q} | {cap - q} |")
        verdict[grp] = {"positions": len(pg), "cap": cap, "quota": q}
    return lines, verdict


def breakdowns(snv: list[dict], indel: list[dict]) -> list[str]:
    lines = []
    elig = [r for r in snv if r["eligible"] == 1]
    nc = [r for r in elig if r["group"] == "noncoding"]
    lines += ["## SNV non-coding group by CADD class", "", "| CADD class | Eligible variants | Positions |", "|---|---|---|"]
    for c in ("splice_region", "synonymous", "utr", "intronic"):
        rows = [r for r in nc if r["class"] == c]
        lines.append(f"| {c} | {len(rows)} | {len(positions_by_gene(rows))} |")
    mois = ["mono", "biallelic", "xl_mono", "xl_bi", "unknown", "none"]
    for name, pool, groups in (("SNV", snv, Q200), ("Indel", indel, Q120)):
        e = [r for r in pool if r["eligible"] == 1]
        lines += ["", f"## {name} groups by inheritance mode (eligible positions)", "",
                  "| Group | " + " | ".join(mois) + " |", "|---|" + "---|" * len(mois)]
        for grp in groups:
            cells = []
            for m in mois:
                cells.append(str(len(positions_by_gene([r for r in e if r["group"] == grp and r["moi"] == m]))))
            lines.append(f"| {grp} | " + " | ".join(cells) + " |")
    lines += ["", "## ClinVar consequence against the SNV class (all SNVs)", "",
              "| ClinVar consequence | " + " | ".join(SNV_CLASS_ORDER) + " |", "|---|" + "---|" * len(SNV_CLASS_ORDER)]
    xt = Counter((r["consequence"], r["class"]) for r in snv)
    for cv in sorted({r["consequence"] for r in snv}, key=lambda c: -sum(n for (a, _), n in xt.items() if a == c)):
        lines.append(f"| {cv} | " + " | ".join(str(xt[(cv, c)]) for c in SNV_CLASS_ORDER) + " |")
    lines += ["", "## ClinVar consequence against the indel class (all indels)", "",
              "| ClinVar consequence | " + " | ".join(INDEL_CLASS_ORDER) + " |",
              "|---|" + "---|" * len(INDEL_CLASS_ORDER)]
    xt = Counter((r["consequence"], r["class"]) for r in indel)
    for cv in sorted({r["consequence"] for r in indel}, key=lambda c: -sum(n for (a, _), n in xt.items() if a == c)):
        lines.append(f"| {cv} | " + " | ".join(str(xt[(cv, c)]) for c in INDEL_CLASS_ORDER) + " |")
    return lines


SNV_CLASS_ORDER = ["nonsense", "missense", "canonical_splice", "splice_region", "synonymous", "utr", "intronic", "other"]
INDEL_CLASS_ORDER = ["frameshift", "inframe", "canonical_splice", "utr", "intronic", "other"]


def summary(sv: dict, iv: dict) -> list[str]:
    lines = ["## Feasibility", ""]
    for grp, v in sv.items():
        q2, q5 = Q200[grp], Q500[grp]
        lo, med, _ = v["sim"]
        sure = v["new"] >= q2 + q5 and v["cap500"] - min(q2, v["cap500"]) >= q5
        if lo >= 0:
            word = "feasible"
        elif med >= 0:
            word = "feasible in most draws, short in some"
        else:
            word = "not feasible as written"
        lines.append(f"- SNV {grp}: {word}. {v['positions']} positions for {q2 + q5} together "
                     f"(spare {v['positions'] - q2 - q5}); {v['new']} not first-run for the 500's {q5} "
                     f"(spare {v['new'] - q5}); under the cap the 500 can draw {v['cap500']} "
                     f"(spare {v['cap500'] - q5}). After a simulated 200, the spare for the 500 is {lo} at worst "
                     f"and {med} at the median. "
                     + ("Feasible even if the 200 took only non-first-run positions." if sure else
                        "Not guaranteed if the 200 took only non-first-run positions.")
                     + (f" The cap bites: {v['over']} gene{'s' if v['over'] > 1 else ''} with more than 10 positions "
                        f"({', '.join(f'{g} {n}' for g, n in v['top'] if n > GENE_CAP)})." if v["over"] else
                        " The cap does not bite."))
    for grp, v in iv.items():
        word = "feasible" if v["cap"] >= v["quota"] else "not feasible as written"
        lines.append(f"- Indel {grp}: {word}. {v['positions']} positions, {v['cap']} under the cap, for "
                     f"{v['quota']} (spare {v['cap'] - v['quota']}).")
    return lines + [""]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cadd-workdir", type=Path, default=Path(tempfile.gettempdir()) / "rdv2_pool_cadd",
                    help="scratch folder for CADD region files and the downloaded tabix index")
    a = ap.parse_args()
    POOL.mkdir(parents=True, exist_ok=True)
    snv_rows, indel_rows = read_tsv(SNV_IN), read_tsv(INDEL_IN)
    snv_rows = [r for r in snv_rows if common.is_snv(r["ref"], r["alt"])]
    indel_rows = [r for r in indel_rows if not common.is_snv(r["ref"], r["alt"])]
    log(f"pools: {len(snv_rows)} SNVs, {len(indel_rows)} indels")
    snv_keys, indel_keys = [key_of(r) for r in snv_rows], [key_of(r) for r in indel_rows]

    cadd_rows = cadd(snv_keys, a.cadd_workdir)
    missing = [k for k in snv_keys if k not in cadd_rows]
    no_mane = sum(1 for k in snv_keys if k in cadd_rows and cadd_rows[k]["mane"] != "True")
    log(f"CADD: {len(missing)} SNVs with no row, {no_mane} with no MANE row (most severe over all rows used)")
    af = gnomad(snv_keys + indel_keys)
    found = Counter((common.is_snv(k[2], k[3]), af[k]["found"]) for k in snv_keys + indel_keys)
    log(f"gnomAD found: SNV {found[(True, 'True')]} of {len(snv_keys)}, indel {found[(False, 'True')]} of "
        f"{len(indel_keys)}")

    snv_class = {k: SNV_CLASS.get(cadd_rows.get(k, {}).get("cadd_consequence", ""), "other") for k in snv_keys}
    indel_class = {k: INDEL_CLASS[c] for k, c in annotate_indels.annotate(INDEL_IN, INDEL_CLASSES).items()}
    # Section 5: no variant with REF and ALT both longer than 1 bp in the indel pool.
    multi = [r for r in indel_rows if len(r["ref"]) > 1 and len(r["alt"]) > 1]
    multi_by_class = Counter(indel_class[key_of(r)] for r in multi)
    indel_rows = [r for r in indel_rows if not (len(r["ref"]) > 1 and len(r["alt"]) > 1)]
    log(f"indel pool: {len(multi)} variants with REF and ALT both longer than 1 bp left out "
        f"({', '.join(f'{c} {n}' for c, n in multi_by_class.most_common())}); {len(indel_rows)} indels left")
    index = common.GeneIndex()
    moi = {r["gene"]: r["moi"] for r in read_tsv(common.DATA / "gene_moi.tsv")}
    snv = build("snv", snv_rows, snv_class, cadd_rows, af, index, moi)
    indel = build("indel", indel_rows, indel_class, {}, af, index, moi)
    for name, pool in (("snv_pool.tsv", snv), ("indel_pool.tsv", indel)):
        write_tsv(POOL / name, list(pool[0]), ([r[c] for c in pool[0]] for r in pool))
        log(f"wrote {POOL / name}: {len(pool)} rows, {sum(r['eligible'] for r in pool)} eligible")

    first = first_run_positions()
    s_lines, sv = snv_section(snv, first)
    i_lines, iv = indel_section(indel, multi_by_class)
    expect = {"missense_variant": "missense", "nonsense": "nonsense", "splice_donor_variant": "canonical_splice",
              "splice_acceptor_variant": "canonical_splice"}
    differ = Counter(r["consequence"] for r in snv if r["consequence"] in expect and r["class"] != expect[r["consequence"]])
    unknown = Counter(r["gene"] for r in snv + indel if r["moi"] == "unknown" and r["eligible"] == 1)
    notes = ["## Notes", "",
             "- CADD class against ClinVar's consequence on the whole SNV pool: "
             + ", ".join(f"{n} of {sum(1 for r in snv if r['consequence'] == c)} {c}"
                         for c, n in sorted(differ.items())) + " differ (table below). PROTOCOL.md section 4 "
             "reports 0 on its smaller check.",
             f"- Inheritance mode unknown (PanelApp gives no mode): {', '.join(sorted(unknown))}, "
             f"{sum(unknown.values())} eligible variants. Section 3 defines no diagnostic genotype for these.", ""]
    md = ["# Diagnostic pool counts", "",
          "Made by `pool.py` (PROTOCOL.md sections 3 to 5). Counts are of eligible variants (full gnomAD v3.1.2 "
          "genome grpmax below 0.1%) and of their distinct positions, since a diagnostic position is used once.",
          f"SNVs with no CADD row: {len(missing)}; with no MANE row in CADD (class from all rows): {no_mane}.", "",
          *summary(sv, iv), *notes, *s_lines, "", *i_lines, "", *breakdowns(snv, indel), ""]
    COUNTS_MD.parent.mkdir(parents=True, exist_ok=True)
    COUNTS_MD.write_text("\n".join(md))
    log(f"wrote {COUNTS_MD}")
    print("\n".join(summary(sv, iv)))
    log("DONE")


if __name__ == "__main__":
    main()
