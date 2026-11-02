"""People of the second run: sex, family, gnomAD release status, and the 50 new replication people.

Implements PROTOCOL.md section 1 (family and gnomAD release), section 3 (sex) and section 5 (replication people).

- Sex, father and mother come from the 1000 Genomes pedigree (1kGP.3202_samples.pedigree_info.txt; 1 male, 2 female).
- A person's family is the person, their parents, their children and their full siblings (same father and same
  mother, both known). The `family` column lists the gnomAD IDs of the members that gnomAD released.
- gnomAD release status and group (gnomad_population_inference.pop) come from the HGDP + 1KG subset sample metadata.
  gnomAD IDs are matched after stripping a trailing A or B (NA12546B, NA12830A, NA18874A).
- Main set: the first run's 50 (rare-disease/data/samples50.txt).
- Replication set: the 2,504 unrelated samples (the 3,202 minus the 698 in the related sequence index), minus both
  first-run sets of 50, minus anyone sharing a 20130606 family ID with them or linked to them in the pedigree as
  parent, child or full sibling. The rest are sorted, shuffled with random.Random(SEEDS["rep_people"]) and the first
  50 taken. Written sorted to data/rep/samples50.txt.

Downloads go to data/gnomad_subset/ (a local cache). Writes data/people.tsv and data/rep/samples50.txt.

  uv run -q --python 3.12 --with-requirements requirements.txt python people.py

Derived from rare-disease/prepare_data.sh (the replication draw) at amp-2026-v4.
"""
from __future__ import annotations

import csv
import gzip
import json
import random
import re
import subprocess
from collections import Counter, defaultdict

from common import DATA, SEEDS, V1_DATA

CACHE = DATA / "gnomad_subset"
KGP = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
URLS = {
    "pedigree": f"{KGP}/working/1kGP.3202_samples.pedigree_info.txt",
    "ped": f"{KGP}/20130606_g1k_3202_samples_ped_population.txt",
    "related": f"{KGP}/1000G_698_related_high_coverage.sequence.index",
    "meta": "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/"
            "gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz",
}


def fetch(name: str):
    path = CACHE / URLS[name].rsplit("/", 1)[1]
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-sSfL", "--retry", "5", "-o", str(path), URLS[name]], check=True)
    return path


def gnomad_meta() -> dict[str, dict]:
    """1000 Genomes ID (A or B suffix stripped) to gnomAD ID, release status and inferred group."""
    csv.field_size_limit(10**9)
    out = {}
    with gzip.open(fetch("meta"), "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if json.loads(r["hgdp_tgp_meta"])["project"] != "1000 Genomes":
                continue
            inferred = r["gnomad_population_inference"]
            pop = json.loads(inferred).get("pop") or "" if inferred.startswith("{") else ""  # "NA" if not inferred
            sample = re.sub(r"[AB]$", "", r["s"])
            assert sample not in out, sample
            out[sample] = {"gnomad_id": r["s"], "gnomad_release": r["gnomad_release"] == "true", "gnomad_pop": pop}
    return out


def pedigree() -> dict[str, dict]:
    with fetch("pedigree").open() as h:
        rows = [line.split() for line in h][1:]
    return {s: {"father": f, "mother": m, "sex": {"1": "male", "2": "female"}[x]} for s, f, m, x in rows}


def family(sample: str, ped: dict[str, dict]) -> set[str]:
    """The person, their parents, children and full siblings (pedigree file)."""
    me = ped[sample]
    fam = {sample} | {p for p in (me["father"], me["mother"]) if p != "0"}
    for s, r in ped.items():
        if sample in (r["father"], r["mother"]):
            fam.add(s)
        if me["father"] != "0" and me["mother"] != "0" and (r["father"], r["mother"]) == (me["father"], me["mother"]):
            fam.add(s)
    return fam


def main() -> None:
    ped = pedigree()
    meta = gnomad_meta()
    with fetch("ped").open() as h:
        fam_id = {r[1]: r[0] for r in (line.split() for line in list(h)[1:])}
    with fetch("related").open() as h:
        related = {line.split("\t")[9] for line in h if not line.startswith("#")}
    main_set = (V1_DATA / "samples50.txt").read_text().split()
    v1_rep = (V1_DATA / "replication" / "samples50.txt").read_text().split()
    assert len(ped) == 3202 and set(ped) == set(fam_id) and set(meta) == set(ped), "sample lists differ"

    # Replication people (PROTOCOL.md section 5).
    unrelated = sorted(set(ped) - related)
    used = set(main_set) | set(v1_rep)
    excluded = set(used)
    used_fam_ids = {fam_id[s] for s in used}
    excluded |= {s for s in ped if fam_id[s] in used_fam_ids}
    for s in used:
        excluded |= family(s, ped)
    candidates = [s for s in unrelated if s not in excluded]
    rng = random.Random(SEEDS["rep_people"])
    rng.shuffle(candidates)
    rep = sorted(candidates[:50])
    (DATA / "rep").mkdir(parents=True, exist_ok=True)
    (DATA / "rep" / "samples50.txt").write_text("\n".join(rep) + "\n")

    with (DATA / "people.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "set", "sex", "gnomad_id", "gnomad_release", "gnomad_pop", "family"])
        for set_name, samples in (("main", sorted(main_set)), ("rep", rep)):
            for s in samples:
                m = meta[s]
                fam = sorted(meta[x]["gnomad_id"] for x in family(s, ped) if meta[x]["gnomad_release"])
                w.writerow([s, set_name, ped[s]["sex"], m["gnomad_id"], str(m["gnomad_release"]).lower(),
                            m["gnomad_pop"], ",".join(fam)])

    # Report.
    print(f"unrelated samples: {len(unrelated)}; left out (first-run sets, family IDs, pedigree links): "
          f"{len(set(unrelated) & excluded)}; candidates: {len(candidates)}")
    print(f"replication people overlapping the first run: {len(set(rep) & used)}")
    for set_name, samples in (("main", main_set), ("rep", rep)):
        rel = [s for s in samples if meta[s]["gnomad_release"]]
        fam_sizes = Counter(sum(meta[x]["gnomad_release"] for x in family(s, ped)) for s in samples)
        print(f"{set_name}: {len(samples)} people, {Counter(ped[s]['sex'] for s in samples)}, "
              f"released in gnomAD {len(rel)}, groups {dict(Counter(meta[s]['gnomad_pop'] for s in rel))}, "
              f"released family members per person {dict(sorted(fam_sizes.items()))}")
    kgp_rel = defaultdict(int)
    for s, m in meta.items():
        if m["gnomad_release"]:
            kgp_rel[m["gnomad_pop"]] += 1
    print(f"released 1000 Genomes samples: {sum(kgp_rel.values())}, by group {dict(sorted(kgp_rel.items()))}")


if __name__ == "__main__":
    main()
