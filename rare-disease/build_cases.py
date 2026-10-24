"""Build rare-disease benchmark cases: a real 1000 Genomes background plus one ClinVar diagnostic variant.

Each case is one 1000 Genomes individual's rare SNVs (1000G AF below 1%) inside the intellectual-disability
panel regions, with one ClinVar pathogenic or likely pathogenic SNV (2+ stars) added as the diagnostic variant.
Background variants that are themselves ClinVar P/LP are dropped so that each case has one correct answer.

Studies:
  snv     SNV diagnostics, SNV background (outputs in data/)
  indel   indel diagnostics (up to 50 bp), SNV and indel background (outputs in data/indel/)
  cancer  PanelApp 245 Adult solid tumours cancer susceptibility: SNV and indel diagnostics and background
          (inputs and outputs in data/cancer/)
  replication  500 new SNV cases to test the snv results on held-out data: diagnostics not at any position used
          as an snv diagnostic, the snv consequence mix scaled to 500, 50 new unrelated 1000 Genomes backgrounds
          (inputs and outputs in data/replication/)

Limitations, stated on the slides: the diagnostic variant is inserted, not observed in that person; inheritance
and phenotype are ignored; ClinVar labels overlap with what some scores were trained on.

  python rare-disease/build_cases.py [snv|indel|cancer|replication]
"""
import csv
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SEED = 20261016
MAX_AF = 0.01
MAX_INDEL = 51  # 50 bp change plus the VCF anchor base
# Diagnostic variants per consequence group. SNV: roughly what reaches a diagnosis in rare disease. Indel: frameshifts
# dominate ClinVar, so the rarer groups are oversampled to give each group enough cases.
STUDIES = {
    "snv": {"clinvar": "clinvar_plp2.tsv", "out": ".", "seed": SEED, "types": 'TYPE="snp"',
            "quota": {"missense": 60, "canonical_splice": 40, "nonsense": 30, "noncoding_or_synonymous": 70}},
    "indel": {"clinvar": "clinvar_indel_plp2.tsv", "out": "indel", "seed": SEED + 1,
              "types": '(TYPE="snp" || TYPE="indel")',
              "quota": {"frameshift": 40, "inframe": 30, "canonical_splice": 30, "noncoding": 20}},
    "cancer": {"clinvar": "cancer/clinvar_all_plp2.tsv", "out": "cancer", "kgp": "cancer/kgp", "seed": SEED + 2,
               "types": '(TYPE="snp" || TYPE="indel")',
               "quota": {"missense": 40, "nonsense": 25, "frameshift": 25, "canonical_splice": 25, "noncoding": 25}},
    "replication": {"clinvar": "clinvar_plp2.tsv", "out": "replication", "kgp": "replication/kgp", "seed": SEED + 7,
                    "types": 'TYPE="snp"', "exclude": "cases.tsv",
                    "quota": {"missense": 150, "canonical_splice": 100, "nonsense": 75, "noncoding_or_synonymous": 175}},
}
# Most severe ClinVar consequence across transcripts. non-coding_transcript_variant is left out: in these genes it
# usually names a non-coding isoform.
NONCODING = {"intron_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "synonymous_variant"}


def group(consequence: str, study: str) -> str | None:
    if study == "cancer":
        if consequence in ("missense_variant", "nonsense", "frameshift_variant"):
            return {"missense_variant": "missense", "nonsense": "nonsense", "frameshift_variant": "frameshift"}[consequence]
        if consequence in ("splice_donor_variant", "splice_acceptor_variant"):
            return "canonical_splice"
        return "noncoding" if consequence in NONCODING else None
    if study == "indel":
        if consequence == "frameshift_variant":
            return "frameshift"
        if consequence in ("inframe_deletion", "inframe_insertion", "inframe_indel"):
            return "inframe"
        if consequence in ("splice_donor_variant", "splice_acceptor_variant"):
            return "canonical_splice"
        if consequence in ("intron_variant", "5_prime_UTR_variant", "3_prime_UTR_variant"):
            return "noncoding"
        return None
    if consequence == "missense_variant":
        return "missense"
    if consequence in ("splice_donor_variant", "splice_acceptor_variant"):
        return "canonical_splice"
    if consequence == "nonsense":
        return "nonsense"
    if consequence in NONCODING:
        return "noncoding_or_synonymous"
    return None


def background(types: str, kgp: Path) -> dict[str, list[tuple]]:
    """Rare variants carried by each sample: sample -> [(chrom, pos, ref, alt, af)]."""
    per_sample = defaultdict(list)
    for bcf in sorted(kgp.glob("*.bcf")):
        samples = subprocess.run(["bcftools", "query", "-l", str(bcf)], capture_output=True, text=True,
                                 check=True).stdout.split()
        out = subprocess.run(["bcftools", "query", "-i", f'{types} && AF<{MAX_AF}',
                              "-f", "%CHROM\t%POS\t%REF\t%ALT\t%AF[\t%GT]\n", str(bcf)],
                             capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            chrom, pos, ref, alt, af, *gts = line.split("\t")
            if "," in alt or max(len(ref), len(alt)) > MAX_INDEL or not set(ref + alt) <= set("ACGT"):
                continue
            for s, gt in zip(samples, gts):
                if "1" in gt:
                    per_sample[s].append((chrom, int(pos), ref, alt, float(af)))
    return per_sample


def main(study: str) -> None:
    cfg = STUDIES[study]
    out = DATA / cfg["out"]
    out.mkdir(exist_ok=True)
    rng = random.Random(cfg["seed"])
    with (DATA / cfg["clinvar"]).open() as h:
        clinvar = list(csv.DictReader(h, delimiter="\t"))
    plp_keys = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in clinvar}

    # Held-out diagnostics: no position used as a diagnostic in the study named by "exclude".
    used = set()
    if "exclude" in cfg:
        with (DATA / cfg["exclude"]).open() as h:
            used = {(r["chrom"], int(r["pos"])) for r in csv.DictReader(h, delimiter="\t") if r["diagnostic"] == "1"}
    by_group = defaultdict(list)
    for r in clinvar:
        g = group(r["consequence"], study)
        if g and (r["chrom"], int(r["pos"])) not in used:
            by_group[g].append(r)
    diagnostic = []
    for g, n in cfg["quota"].items():
        for r in rng.sample(by_group[g], n):
            diagnostic.append({**r, "group": g})
    rng.shuffle(diagnostic)

    per_sample = background(cfg["types"], DATA / cfg.get("kgp", "kgp"))
    samples = sorted(per_sample)
    rows = []
    for i, d in enumerate(diagnostic):
        case = f"case{i + 1:03d}"
        sample = samples[i % len(samples)]
        rows.append([case, sample, d["chrom"], d["pos"], d["ref"], d["alt"], 1, d["group"], d["gene"],
                     d["clinvar_id"]])
        seen = {(d["chrom"], int(d["pos"]), d["ref"], d["alt"])}
        for chrom, pos, ref, alt, af in per_sample[sample]:
            key = (chrom, pos, ref, alt)
            if key in plp_keys or key in seen:
                continue
            seen.add(key)
            rows.append([case, sample, chrom, pos, ref, alt, 0, "", "", ""])

    with (out / "cases.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["case", "sample", "chrom", "pos", "ref", "alt", "diagnostic", "group", "gene", "clinvar_id"])
        w.writerows(rows)
    variants = sorted({(r[2], int(r[3]), r[4], r[5]) for r in rows}, key=lambda v: (v[0], v[1], v[2], v[3]))
    with (out / "case_variants.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt"])
        w.writerows(variants)

    sizes = Counter(r[0] for r in rows)
    print(f"{len(diagnostic)} cases, {len(samples)} samples, {len(variants)} unique variants")
    print(f"variants per case: min {min(sizes.values())}, median {sorted(sizes.values())[len(sizes) // 2]},"
          f" max {max(sizes.values())}")
    print("diagnostic groups:", dict(Counter(d["group"] for d in diagnostic)))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "snv")
