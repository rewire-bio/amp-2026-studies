"""Build the second run's rare-disease cases: a real 1000 Genomes background plus one ClinVar diagnostic variant.

Implements PROTOCOL.md sections 2 (frequency filter), 3 (diagnostic genotype from the mode of inheritance) and 5
(cases, quotas, gene cap, seeds).

Diagnostics come from the classified pools written by pool.py (data/pool/snv_pool.tsv, indel_pool.tsv): eligible
rows (full gnomAD grpmax below 0.1%) in the study's groups. For each group, in the protocol's order, the candidates
are sorted by key and put in a random order with the study's seed (rng.sample over the whole list, since one
rng.sample of the quota cannot respect the cap). They are then taken in that order, skipping a position already used
by this study and a gene that already has 10 diagnostics in the group (the pool's gene), until the quota is met. The
replication study also skips every position used as a diagnostic in this run's SNV cases or in the first run's 200
and 500. A group that cannot fill its quota is reduced to what it has, and the reduction is printed and saved. The
drawn diagnostics are shuffled and assigned round-robin to the set's 50 people, sorted, as in the first run.

Diagnostic zygosity (section 3): hemi for a male on chrX outside the PARs; else hom in biallelic genes and, for
females, in xl_bi genes; else het (mono, xl_mono and unknown, which is GNAS only).

Background: the person's variants (person_variants.tsv, SNVs and indels in every study) with family-left-out grpmax
below 1%, minus ClinVar P/LP variants of any review status (from clinvar.vcf.gz) and the diagnostic itself. Every
study keeps both types because the inheritance filter (section 3) counts the person's variants of both types in a
gene; case_benchmark.py then ranks only variants of the study's type (SNVs for snv and replication, indels for
indel). The main 0.1% filter is applied at analysis time by case_benchmark.py. The draws do not use the
backgrounds, so keeping the other type does not change them.

Columns of cases.tsv: case, sample, chrom, pos, ref, alt, diagnostic, group, class, gene, moi, zygosity, af,
grpmax, grpmax_no1kg, clinvar_id. gene lists every MANE region gene at the variant (common.GeneIndex), comma-joined
and sorted, except that a diagnostic's pool gene comes first; moi gives each gene's mode in the same order.
Background rows take af, grpmax and grpmax_no1kg from person_variant_af.tsv (af_fam, grpmax_fam, grpmax_no1kg);
diagnostics take af and grpmax from full gnomAD (data/pool/gnomad_full.tsv) and grpmax_no1kg equal to grpmax, as
the person does not carry them. Background class is left empty: it is assigned from CADD at analysis time.
case_variants.tsv holds the distinct keys across the study's cases.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python build_cases.py \\
      --study snv|indel|replication [--cap AF]

--cap replaces the 1% background grpmax cap and adds a kgp_af column (the background variant's 1000 Genomes AF,
empty for diagnostics), which the frequency sweep needs (PROTOCOL.md section 13). The sweep runs it with --cap 0.05
and RDV2_DATA=data/af5, RDV2_RESULTS=results/af5, so the v6 files are not touched (AF_SWEEP_RUN.md).

Writes data/<study>/cases.tsv, data/<study>/case_variants.tsv and results/build_<study>.md.

Derived from rare-disease/build_cases.py at amp-2026-v4.
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import statistics
import subprocess
from collections import Counter, defaultdict

import common

QUOTAS = {
    "snv": {"missense": 60, "canonical_splice": 40, "nonsense": 30, "noncoding": 70},
    "indel": {"frameshift": 40, "inframe": 30, "canonical_splice": 30, "noncoding": 20},
    "replication": {"missense": 150, "canonical_splice": 100, "nonsense": 75, "noncoding": 175},
}
POOL = {"snv": "snv_pool.tsv", "indel": "indel_pool.tsv", "replication": "snv_pool.tsv"}
GENE_CAP = 10
FIRST_RUN_CASES = [common.V1_DATA / "cases.tsv", common.V1_DATA / "replication" / "cases.tsv"]
CLINVAR_VCF = common.V1_DATA / "clinvar.vcf.gz"
PLP_SIG = ("Pathogenic", "Likely_pathogenic")
HEADER = ["case", "sample", *common.KEY, "diagnostic", "group", "class", "gene", "moi", "zygosity", "af", "grpmax",
          "grpmax_no1kg", "clinvar_id"]


def clinvar_plp_any() -> set[tuple]:
    """Keys of every ClinVar P/LP variant, any review status (protocol section 5: background P/LP removed).

    A variant counts when the first part of CLNSIG (before any '|') starts with Pathogenic or Likely_pathogenic, so
    Pathogenic/Likely_pathogenic and the low-penetrance labels count and conflicting classifications do not.
    """
    out = subprocess.run(["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CLNSIG\\n", str(CLINVAR_VCF)],
                         capture_output=True, text=True, check=True).stdout
    keys = set()
    for line in out.splitlines():
        chrom, pos, ref, alt, sig = line.split("\t")
        if sig.split("|")[0].startswith(PLP_SIG):
            keys.add(("chr" + chrom, int(pos), ref, alt))
    return keys


def read_tsv(path) -> list[dict]:
    with open(path, newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def key_of(r: dict) -> tuple:
    return (r["chrom"], int(r["pos"]), r["ref"], r["alt"])


def variant_end(pos: int, ref: str, alt: str) -> int:
    """Last reference base touched, as pool.py uses it (an insertion touches the base after the anchor)."""
    return pos + max(len(ref) - 1, 1 if len(alt) > 1 else 0)


def diagnostic_positions(path) -> set[tuple]:
    return {(r["chrom"], int(r["pos"])) for r in read_tsv(path) if r["diagnostic"] == "1"}


def draw(study: str, pool: list[dict], excluded: set[tuple], rng: random.Random, log: list[str],
         capacity: dict[str, int]) -> list[dict]:
    """Diagnostics per group under the gene cap, one per position; reduce a quota that cannot be filled.

    capacity gets, per group, the positions it could still supply under the cap when its turn came.
    """
    used: set[tuple] = set()
    chosen = []
    for grp, quota in QUOTAS[study].items():
        cands = sorted((r for r in pool if r["eligible"] == "1" and r["group"] == grp
                        and (r["chrom"], int(r["pos"])) not in excluded), key=key_of)
        free = {(r["chrom"], int(r["pos"])): r["gene"] for r in cands if (r["chrom"], int(r["pos"])) not in used}
        capacity[grp] = sum(min(n, GENE_CAP) for n in Counter(free.values()).values())
        order = rng.sample(cands, len(cands))
        per_gene: Counter = Counter()
        taken = []
        for r in order:
            if len(taken) == quota:
                break
            p = (r["chrom"], int(r["pos"]))
            if p in used or per_gene[r["gene"]] >= GENE_CAP:
                continue
            used.add(p)
            per_gene[r["gene"]] += 1
            taken.append({**r, "group": grp})
        if len(taken) < quota:
            log.append(f"Quota reduced: {grp} {quota} to {len(taken)} ({len(cands)} candidate rows).")
        chosen += taken
    rng.shuffle(chosen)
    return chosen


def zygosity(moi: str, sex: str, chrom: str, pos: int) -> str:
    if common.haploid_region(chrom, pos, sex):
        return "hemi"
    if moi == "biallelic" or (moi == "xl_bi" and sex == "female"):
        return "hom"
    return "het"


def summary_stats(values: list[int]) -> str:
    return f"median {statistics.median(values):g} ({min(values)} to {max(values)})"


def main(study: str, cap: float = common.LOOSE_AF) -> None:
    rng = random.Random(common.SEEDS[study])
    people_set = common.STUDIES[study]
    out_dir = common.study_dir(study)
    out_dir.mkdir(parents=True, exist_ok=True)
    log: list[str] = []

    people = {r["sample"]: r for r in read_tsv(common.DATA / "people.tsv") if r["set"] == people_set}
    samples = sorted(people)
    moi_of = {r["gene"]: r["moi"] for r in read_tsv(common.DATA / "gene_moi.tsv")}
    index = common.GeneIndex()
    pool = read_tsv(common.DATA / "pool" / POOL[study])
    full = {key_of(r): r for r in read_tsv(common.DATA / "pool" / "gnomad_full.tsv")}
    plp = clinvar_plp_any()

    excluded: set[tuple] = set()
    if study == "replication":
        for f in [common.study_dir("snv") / "cases.tsv", *FIRST_RUN_CASES]:
            pos = diagnostic_positions(f)
            log.append(f"Excluded {len(pos)} diagnostic positions of {os.path.relpath(f, common.HERE.parent)}.")
            excluded |= pos
        log.append(f"Excluded positions in all: {len(excluded)}.")
    capacity: dict[str, int] = {}
    diagnostics = draw(study, pool, excluded, rng, log, capacity)

    missing_moi: Counter = Counter()

    def genes_at(chrom: str, pos: int, ref: str, alt: str, first: str = "") -> tuple[str, str]:
        names = index.genes(chrom, pos, variant_end(pos, ref, alt))
        genes = sorted({g for n in names for g in n.split(",")})
        if first in genes:
            genes = [first] + [g for g in genes if g != first]
        mois = []
        for g in genes:
            if g not in moi_of:
                missing_moi[g] += 1
            mois.append(moi_of.get(g, "none"))
        return ",".join(genes), ",".join(mois)

    # Backgrounds: person variants joined to their family-left-out frequencies.
    afs = {(r["sample"], *key_of(r)): r for r in read_tsv(common.set_dir(people_set) / "person_variant_af.tsv")}
    per_sample: dict[str, list[dict]] = defaultdict(list)
    for r in read_tsv(common.set_dir(people_set) / "person_variants.tsv"):
        if r["sample"] in people:
            per_sample[r["sample"]].append(r)

    rows = []
    stats = {"bg_before": [], "bg_after": [], "plp_removed": 0, "same_pos": 0, "overlap": 0, "unjoined": 0}
    for i, d in enumerate(diagnostics):
        case = f"case{i + 1:03d}"
        sample = samples[i % len(samples)]
        dkey = key_of(d)
        sex = people[sample]["sex"]
        g = full[dkey]
        gene, moi = genes_at(*dkey, first=d["gene"])
        zyg = zygosity(d["moi"], sex, dkey[0], dkey[1])
        rows.append({"case": case, "sample": sample, "chrom": dkey[0], "pos": dkey[1], "ref": dkey[2],
                     "alt": dkey[3], "diagnostic": 1, "group": d["group"], "class": d["class"], "gene": gene,
                     "moi": moi, "zygosity": zyg, "af": g["af"], "grpmax": d["grpmax"], "grpmax_no1kg": d["grpmax"],
                     "clinvar_id": d["clinvar_id"], "_full_grpmax": g["grpmax"],
                     "_dmoi": d["moi"], "_sex": sex, "_pgene": d["gene"]})
        seen = {dkey}
        before = after = 0
        for v in per_sample[sample]:
            key = key_of(v)
            if key in seen:
                continue
            a = afs.get((sample, *key))
            if a is None:
                stats["unjoined"] += 1
                continue
            if float(a["grpmax_fam"]) >= cap:
                continue
            if key in plp:
                stats["plp_removed"] += 1
                continue
            seen.add(key)
            if key[:2] == dkey[:2]:
                stats["same_pos"] += 1
            if key[0] == dkey[0] and key[1] <= variant_end(*dkey[1:]) and dkey[1] <= variant_end(*key[1:]):
                stats["overlap"] += 1
            gene, moi = genes_at(*key)
            rows.append({"case": case, "sample": sample, "chrom": key[0], "pos": key[1], "ref": key[2],
                         "alt": key[3], "diagnostic": 0, "group": "", "class": "", "gene": gene, "moi": moi,
                         "zygosity": v["zygosity"], "af": a["af_fam"], "grpmax": a["grpmax_fam"],
                         "grpmax_no1kg": a["grpmax_no1kg"], "clinvar_id": "", "kgp_af": v["kgp_af"]})
            before += 1
            after += float(a["grpmax_fam"]) < common.MAIN_AF
        stats["bg_before"].append(before)
        stats["bg_after"].append(after)

    with (out_dir / "cases.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, HEADER + (["kgp_af"] if cap != common.LOOSE_AF else []), delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    variants = sorted({(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows})
    with (out_dir / "case_variants.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(common.KEY)
        w.writerows(variants)

    report(study, rows, samples, stats, missing_moi, log, len(variants), capacity, cap)


def report(study, rows, samples, stats, missing_moi, log, n_variants, capacity, cap=common.LOOSE_AF) -> None:
    diag = [r for r in rows if r["diagnostic"] == 1]
    bg = [r for r in rows if r["diagnostic"] == 0]
    per_person = Counter(r["sample"] for r in diag)
    groups = list(QUOTAS[study])
    gene_counts = Counter((r["group"], r["_pgene"]) for r in diag)

    # Checks.
    checks = {
        "every diagnostic grpmax below 0.1%": all(float(r["grpmax"]) < common.MAIN_AF for r in diag),
        f"every background grpmax below {cap * 100:g}%": all(float(r["grpmax"]) < cap for r in bg),
        "every diagnostic's pool grpmax equals full gnomAD grpmax": all(
            float(r["grpmax"]) == float(r["_full_grpmax"]) for r in diag),
        "no gene over 10 in a group": max(gene_counts.values()) <= GENE_CAP,
        "one diagnostic per position": len({(r["chrom"], r["pos"]) for r in diag}) == len(diag),
        "every hemi row is male, chrX, outside the PARs": all(
            r["chrom"] == "chrX" and not common.in_par(r["chrom"], r["pos"])
            for r in rows if r["zygosity"] == "hemi"),
        "every background row joins to a frequency": stats["unjoined"] == 0,
        "every gene has a mode of inheritance": not missing_moi,
    }
    hemi_samples = {r["sample"] for r in rows if r["zygosity"] == "hemi"}
    sexes = {r["sample"]: r["_sex"] for r in diag}
    checks["every hemi row is male, chrX, outside the PARs"] &= all(sexes.get(s, "male") == "male"
                                                                   for s in hemi_samples)
    checks["no background variant has REF and ALT both longer than 1 bp"] = not any(
        len(r["ref"]) > 1 and len(r["alt"]) > 1 for r in bg)
    if study == "indel":
        checks["no diagnostic has REF and ALT both longer than 1 bp"] = not any(
            len(r["ref"]) > 1 and len(r["alt"]) > 1 for r in diag)
    else:
        checks["every diagnostic is an SNV"] = all(common.is_snv(r["ref"], r["alt"]) for r in diag)

    lines = [f"# build_cases.py --study {study}", "",
             f"{len(diag)} cases, {len(samples)} people ({common.STUDIES[study]} set), {len(rows)} case rows, "
             f"{n_variants} distinct variants. Seed {common.SEEDS[study]}.", "",
             "Cases per person: " + ", ".join(f"{n} people with {k}" for k, n in
                                              sorted(Counter(per_person.values()).items())) + ".", "",
             f"Background variants removed as ClinVar P/LP (any review status): {stats['plp_removed']}. Background variants at "
             f"the diagnostic's position: {stats['same_pos']}; overlapping the diagnostic's reference bases (VCF POS to "
             f"the last base touched): {stats['overlap']}.", "",
             "Background per case, SNVs and indels counted apart. Both are kept for the inheritance filter; only the "
             f"study's type ({'indels' if study == 'indel' else 'SNVs'}) is ranked. The inheritance filter is applied "
             "at analysis time and is not counted here.", "",
             f"| Type | Build filter (grpmax below {cap * 100:g}%) | Main filter (grpmax below 0.1%) | Total rows |",
             "|---|---|---|---|"]
    cases = [r["case"] for r in diag]
    for name, want_snv in (("SNVs", True), ("Indels", False)):
        sel = [r for r in bg if common.is_snv(r["ref"], r["alt"]) == want_snv]
        loose = Counter(r["case"] for r in sel)
        main_f = Counter(r["case"] for r in sel if float(r["grpmax"]) < common.MAIN_AF)
        lines.append(f"| {name} | {summary_stats([loose[c] for c in cases])} | "
                     f"{summary_stats([main_f[c] for c in cases])} | {len(sel)} |")
    lines.append(f"| Both | {summary_stats(stats['bg_before'])} | {summary_stats(stats['bg_after'])} | {len(bg)} |")
    lines.append("")
    lines += ["## Quotas", ""] + log + [
              "Every quota filled." if not any(x.startswith("Quota reduced") for x in log) else "", "",
              "Available: positions the group could still supply under the gene cap when it was drawn (after "
              "exclusions and positions taken by earlier groups).", "",
              "Largest gene: the most diagnostics from one gene in the group (the cap is 10 per group).", "",
              "| Group | Quota | Available | Drawn | Genes | Largest gene | mono | biallelic | xl_mono | xl_bi | unknown "
              "| het | hom | hemi |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for grp in groups:
        d = [r for r in diag if r["group"] == grp]
        moi = Counter(r["_dmoi"] for r in d)
        zyg = Counter(r["zygosity"] for r in d)
        genes = Counter(r["_pgene"] for r in d)
        lines.append(f"| {grp} | {QUOTAS[study][grp]} | {capacity[grp]} | {len(d)} | {len(genes)} | "
                     f"{max(genes.values()) if genes else 0} | "
                     + " | ".join(str(moi[m]) for m in ("mono", "biallelic", "xl_mono", "xl_bi", "unknown"))
                     + " | " + " | ".join(str(zyg[z]) for z in ("het", "hom", "hemi")) + " |")
    moi_all = Counter(r["_dmoi"] for r in diag)
    zyg_all = Counter(r["zygosity"] for r in diag)
    lines.append(f"| all | {sum(QUOTAS[study].values())} | | {len(diag)} | {len({r['_pgene'] for r in diag})} | "
                 "| "
                 + " | ".join(str(moi_all[m]) for m in ("mono", "biallelic", "xl_mono", "xl_bi", "unknown"))
                 + " | " + " | ".join(str(zyg_all[z]) for z in ("het", "hom", "hemi")) + " |")
    lines += ["", "Diagnostic classes: " + ", ".join(f"{c} {n}" for c, n in
                                                     Counter(r["class"] for r in diag).most_common()) + ".",
              "Diagnostics in more than one gene: " + str(sum("," in r["gene"] for r in diag)) + ".",
              "Background zygosity: " + ", ".join(f"{z} {n}" for z, n in
                                                  Counter(r["zygosity"] for r in bg).most_common()) + ".",
              "Background variants with no gene: " + str(sum(r["gene"] == "" for r in bg)) + ".", "",
              "## Checks", ""]
    lines += [f"- {name}: {'yes' if ok else 'NO'}" for name, ok in checks.items()]
    if missing_moi:
        lines.append(f"- genes with no mode of inheritance: {dict(missing_moi)}")
    text = "\n".join(lines) + "\n"
    print(text)
    common.RESULTS.mkdir(parents=True, exist_ok=True)
    (common.RESULTS / f"build_{study}.md").write_text(text)
    if not all(checks.values()):
        raise SystemExit("a check failed")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--study", required=True, choices=list(QUOTAS))
    ap.add_argument("--cap", type=float, default=common.LOOSE_AF,
                    help="background grpmax cap (default 1%%); another cap adds the kgp_af column")
    a = ap.parse_args()
    main(a.study, a.cap)
