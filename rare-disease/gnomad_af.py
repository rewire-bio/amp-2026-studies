"""Look up gnomAD genome allele frequencies for every case variant, for the frequency-filtered rerun.

A diagnostic pipeline removes variants that are common in gnomAD before it ranks anything. The case backgrounds are
rare in 1000 Genomes (AF below 1%) but were never filtered on gnomAD, and AlphaGenome AVI is trained to tell
gnomAD-rare from gnomAD-common variants. This script records each variant's gnomAD genome AF and its grpmax AF (the
highest of the non-bottlenecked genetic ancestry groups: afr, amr, eas, mid, nfe, sas) so that
case_benchmark.py --gnomad-max-af can drop what a pipeline would drop. 1000 Genomes samples are part of gnomAD, so
nearly every background variant is found; a variant that is not found gets AF 0 and is kept.

Frequencies come from MyVariant.info's gnomad_genome field (gnomAD v3.1.2 genomes, the same 76,215 genomes as
gnomAD v4.1 genomes; spot checks match the v4.1 VCF AF exactly), queried in batches of 1,000 by HGVS id. Indels use
VCF-style (left-aligned) positions.

  python rare-disease/gnomad_af.py      # writes data/gnomad_af.tsv
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
URL = "https://myvariant.info/v1/variant?assembly=hg38&fields=gnomad_genome.af"
KEY = ["chrom", "pos", "ref", "alt"]
STUDIES = ["", "replication", "cancer", "indel"]
GROUPS = ["afr", "amr", "eas", "mid", "nfe", "sas"]


def hgvs(chrom: str, pos: int, ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}{ref}>{alt}"
    if alt == ref[0] and len(ref) > 1:  # deletion
        s, e = pos + 1, pos + len(ref) - 1
        return f"{chrom}:g.{s}del" if s == e else f"{chrom}:g.{s}_{e}del"
    if ref == alt[0] and len(alt) > 1:  # insertion
        return f"{chrom}:g.{pos}_{pos + 1}ins{alt[1:]}"
    e = pos + len(ref) - 1
    return f"{chrom}:g.{pos}_{e}delins{alt}" if e > pos else f"{chrom}:g.{pos}delins{alt}"


def fetch(ids: list[str]) -> dict[str, tuple[float, float]]:
    body = urllib.parse.urlencode({"ids": ",".join(ids)}).encode()
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(URL, data=body), timeout=120) as r:
                hits = json.load(r)
            break
        except OSError:
            time.sleep(5 * (attempt + 1))
    else:
        raise SystemExit("MyVariant.info did not answer")
    out = {}
    for h in hits:
        af = (h.get("gnomad_genome") or {}).get("af")
        if isinstance(af, list):
            af = af[0]
        if af:
            out[h["query"]] = (float(af.get("af", 0.0)), max(float(af.get(f"af_{g}", 0.0)) for g in GROUPS))
    return out


def main() -> None:
    variants = pd.concat([pd.read_csv(DATA / s / "cases.tsv", sep="\t", usecols=KEY) for s in STUDIES]).drop_duplicates()
    variants["id"] = [hgvs(*r) for r in variants[KEY].itertuples(index=False)]
    ids = variants["id"].tolist()
    found: dict[str, tuple[float, float]] = {}
    for i in range(0, len(ids), 1000):
        found |= fetch(ids[i:i + 1000])
        print(f"{min(i + 1000, len(ids))} of {len(ids)} queried, {len(found)} found", flush=True)
    variants["gnomad_found"] = variants["id"].isin(found)
    variants["gnomad_af"] = [found.get(i, (0.0, 0.0))[0] for i in variants["id"]]
    variants["gnomad_af_grpmax"] = [found.get(i, (0.0, 0.0))[1] for i in variants["id"]]
    variants.to_csv(DATA / "gnomad_af.tsv", sep="\t", index=False)
    snv = (variants["ref"].str.len() == 1) & (variants["alt"].str.len() == 1)
    print(f"{len(variants)} variants; found in gnomAD genomes: SNVs {variants.loc[snv, 'gnomad_found'].mean():.1%}, "
          f"indels {variants.loc[~snv, 'gnomad_found'].mean():.1%}")


if __name__ == "__main__":
    main()
