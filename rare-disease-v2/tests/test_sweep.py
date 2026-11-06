"""The frequency sweep (PROTOCOL.md section 13) on synthetic cases.

Writes a synthetic copy of the data in a scratch folder: 50 main and 50 replication people, an snv study (200 cases)
and a replication study (100 cases), each built twice, as the v6 build (background grpmax and 1000 Genomes AF below
1%) and as the 5%-cap build in data/af5/ (below 5%, with a kgp_af column, the same rows in the same order below 1%).
Scores are random, with AVI lower for common variants as the sweep's hypothesis has it. Then runs case_benchmark.py
--analysis main and bg1pct on the v6 copy and --analysis sweep on the 5% copy (RDV2_DATA, RDV2_RESULTS and RDV2_EVO2
point there) and checks that:

  - check_af5.py finds the two builds identical below the 1% caps, and finds one changed frequency;
  - the 0.1% and 1% sweep rows reproduce the main and bg1pct results (the sweep's own check and one done here);
  - a case is dropped exactly when its diagnostic's grpmax is at or above the threshold;
  - Evo 2 is reported at 1% and below only, and the median list never shrinks as the threshold rises;
  - competitors by grpmax band add up to competitors() for AVI and CADD;
  - results/analysis_sweep.md is written, with every threshold and no em or en dash;
  - figures/slide_charts.py draws af_sweep and af_sweep_gain, PNG and SVG, from the sweep result.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python tests/test_sweep.py

Exit status 0 when every check passes.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import case_benchmark as cb  # noqa: E402
import common  # noqa: E402

KEY = common.KEY
SEED = 20261106
CONSEQUENCES = {"nonsense": "STOP_GAINED", "missense": "NON_SYNONYMOUS", "canonical_splice": "CANONICAL_SPLICE",
                "intronic": "INTRONIC", "splice_region": "SPLICE_SITE"}
GROUP = {"nonsense": "nonsense", "missense": "missense", "canonical_splice": "canonical_splice",
         "intronic": "noncoding", "splice_region": "noncoding"}
GENES = [(f"G{i}", "biallelic" if i % 3 == 0 else "mono") for i in range(30)]


def synthetic(root: Path) -> dict:
    """Write the synthetic data under root/data and root/evo2; return the diagnostics' grpmax per study."""
    rng = np.random.default_rng(SEED)
    data = root / "data"
    people = [(f"M{i:02d}", "main") for i in range(50)] + [(f"R{i:02d}", "rep") for i in range(50)]
    (data / common.SWEEP_DATA).mkdir(parents=True)
    pd.DataFrame({"sample": [p for p, _ in people], "set": [s for _, s in people],
                  "sex": ["female" if i % 2 else "male" for i in range(len(people))]}).to_csv(
        data / "people.tsv", sep="\t", index=False)

    pos = iter(range(900_000_001, 999_999_999, 7))  # beyond chr1's end, so no real key is reused
    scores, gnomad, evo2, diag_af = [], [], [], {}

    def variant(cls: str, af: float, diagnostic: bool) -> dict:
        p = next(pos)
        snv = diagnostic or rng.random() > 0.08
        ref, alt = ("A", "G") if snv else ("AT", "A")
        commonness = np.log10(af + 1e-6)  # -6 (absent) to about -1.3 (5%)
        hit = 2.0 if diagnostic else 0.0
        row = {"chrom": "chr1", "pos": p, "ref": ref, "alt": alt}
        scores.append(row | {
            "cadd_phred": round(float(rng.normal(14 + 6 * hit, 6)), 3),
            "cadd_consequence": CONSEQUENCES.get(cls, cls), "spliceai_max": round(float(rng.beta(1, 8 - 2 * hit)), 3),
            "revel": round(float(rng.beta(2 + hit, 4)), 3) if cls == "missense" else np.nan,
            "alphamissense": round(float(rng.beta(2 + hit, 4)), 3) if cls == "missense" else np.nan,
            "avi": round(float(rng.normal(hit - 0.6 * (commonness + 6) / 4.7, 0.8)), 4),
            "gpn_msa": round(float(rng.normal(hit, 2)), 3), "gpnstar_v": round(float(rng.normal(hit, 2)), 3),
            "gpnstar_m": round(float(rng.normal(hit, 2)), 3)})
        evo2.append(row | {"evo2_delta": round(float(rng.normal(-hit, 1)), 4)})
        return row

    for study, people_set, n_cases, mix in (
            ("snv", "main", 200, {"nonsense": 30, "missense": 60, "canonical_splice": 40, "intronic": 40,
                                  "splice_region": 30}),
            ("replication", "rep", 100, {"nonsense": 15, "missense": 30, "canonical_splice": 20, "intronic": 20,
                                         "splice_region": 15})):
        samples = sorted(p for p, s in people if s == people_set)
        backgrounds = {}
        for s in samples:  # one background per person, shared by their cases, in a fixed order
            rows = []
            for _ in range(int(rng.integers(25, 45))):
                grp = 0.0 if rng.random() < 0.05 else float(10 ** rng.uniform(-5, np.log10(0.0499)))
                kgp = min(0.0499, grp * float(rng.uniform(0.5, 2.5)))
                cls = rng.choice(["intronic", "missense", "SYNONYMOUS", "3PRIME_UTR", "splice_region"],
                                 p=[0.5, 0.2, 0.15, 0.1, 0.05])
                gene, moi = GENES[int(rng.integers(len(GENES)))]
                rows.append(variant(str(cls), grp, False) | {
                    "gene": gene, "moi": moi, "zygosity": "hom" if rng.random() < 0.1 else "het",
                    "af": grp / 2, "grpmax": grp, "grpmax_no1kg": grp, "kgp_af": round(kgp, 6)})
            backgrounds[s] = rows
        classes = [c for c, n in mix.items() for _ in range(n)]
        rng.shuffle(classes)
        rows, afs = [], []
        for i, cls in enumerate(classes[:n_cases]):
            case, sample = f"case{i + 1:03d}", samples[i % len(samples)]
            af = float(rng.choice([0.0, 0.00005, 0.0002, 0.0007], p=[0.8, 0.1, 0.06, 0.04]))
            afs.append(af)
            gene, moi = GENES[int(rng.integers(len(GENES)))]
            d = variant(cls, af, True)
            gnomad.append(d | {"found": int(af > 0 or rng.random() < 0.1), "ac": int(af > 0)})
            rows.append({"case": case, "sample": sample, **d, "diagnostic": 1, "group": GROUP[cls], "class": cls,
                         "gene": gene, "moi": moi, "zygosity": "hom" if moi == "biallelic" else "het",
                         "af": af / 2, "grpmax": af, "grpmax_no1kg": af, "clinvar_id": str(1000 + i), "kgp_af": ""})
            rows += [{"case": case, "sample": sample, **b, "diagnostic": 0, "group": "", "class": "",
                      "clinvar_id": ""} for b in backgrounds[sample]]
        diag_af[study] = afs
        af5 = pd.DataFrame(rows)
        bg = af5["diagnostic"] == 0
        v6 = af5[~bg | ((af5["grpmax"] < common.LOOSE_AF) & (af5["kgp_af"].replace("", np.nan).astype(float)
                                                             < common.KGP_MAX_AF))].drop(columns="kgp_af")
        for folder, cases in ((data / study, v6), (data / common.SWEEP_DATA / study, af5)):
            folder.mkdir(parents=True)
            cases.to_csv(folder / "cases.tsv", sep="\t", index=False)
            keys = cases[KEY].drop_duplicates()
            s = pd.DataFrame(scores).merge(keys, on=KEY)
            s[KEY + ["cadd_phred", "cadd_consequence", "spliceai_max"]].to_csv(folder / "scores_cadd.tsv", sep="\t",
                                                                             index=False)
            for source, col in (("revel", "revel"), ("am", "alphamissense"), ("avi", "avi"), ("gpnmsa", "gpn_msa"),
                                ("gpnstar_v", "gpnstar_v"), ("gpnstar_m", "gpnstar_m")):
                s[KEY + [col]].dropna().to_csv(folder / f"scores_{source}.tsv", sep="\t", index=False)

    pd.DataFrame(gnomad)[KEY + ["found", "ac"]].to_csv(data / "gnomad_full.tsv", sep="\t", index=False)
    # Evo 2, as in the GPU run, covers only variants below 1%: the v6 lists.
    e = pd.DataFrame(evo2)
    v6_keys = pd.concat([pd.read_csv(data / s / "cases.tsv", sep="\t")[KEY] for s in ("snv", "replication")])
    (root / "evo2" / "run1").mkdir(parents=True)
    e.merge(v6_keys.drop_duplicates(), on=KEY).to_csv(root / "evo2" / "run1" / "scores.tsv", sep="\t", index=False)
    return diag_af


def run(root: Path, study: str, analysis: str) -> dict:
    env = os.environ | {"RDV2_DATA": str(root / "data"), "RDV2_RESULTS": str(root / "results"),
                        "RDV2_EVO2": str(root / "evo2")}
    subprocess.run([sys.executable, str(HERE / "case_benchmark.py"), "--study", study, "--analysis", analysis],
                   cwd=HERE, env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return json.loads((root / "results" / f"{study}_{analysis}.json").read_text())


def check_bands() -> list[str]:
    """competitors_by_band adds up to competitors(), a missing diagnostic score counting every background variant."""
    cases = pd.DataFrame({
        "case": ["a"] * 5 + ["b"] * 3, "diagnostic": [1, 0, 0, 0, 0, 1, 0, 0],
        "grpmax": [0.0, 0.0005, 0.005, 0.03, 0.0, 0.0, 0.002, 0.04],
        "AVI": [1.0, 2.0, 1.0, 0.5, np.nan, np.nan, 0.1, 3.0]})
    by = cb.competitors_by_band(cases, "AVI")
    problems = []
    if not by.sum(axis=1).equals(cb.competitors(cases, "AVI").astype(by.sum(axis=1).dtype)):
        problems.append(f"bands do not add up: {by.to_dict()}")
    want = {"a": [1, 1, 0], "b": [0, 1, 1]}  # a: 2.0 (below 0.1%) and 1.0 tied (0.1 to 1%); b: diagnostic unscored
    for case, counts in want.items():
        if list(by.loc[case]) != counts:
            problems.append(f"case {case}: {list(by.loc[case])}, expected {counts}")
    return problems


def leaves(x, path=()) -> dict:
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            out |= leaves(v, path + (k,))
        return out
    return {path: x}


def check_identity(root: Path) -> list[str]:
    """check_af5.py passes on the synthetic builds, then fails once one background frequency in af5 is changed."""
    env = os.environ | {"RDV2_DATA": str(root / "data"), "RDV2_RESULTS": str(root / "results")}
    cmd = [sys.executable, str(HERE / "check_af5.py")]
    problems = []
    if subprocess.run(cmd, cwd=HERE, env=env, stdout=subprocess.DEVNULL).returncode != 0:
        problems.append("check_af5.py: the synthetic builds should be identical: "
                        + (root / "results" / "af5_identity.md").read_text()[-600:])
    path = root / "data" / common.SWEEP_DATA / "replication" / "cases.tsv"
    original = path.read_text()
    cases = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    i = cases.index[(cases["diagnostic"] == "0") & (cases["grpmax"].astype(float) < 0.0005)][0]
    cases.loc[i, "grpmax"] = "0.0123"  # now above the 1% cap: the v6 row is missing from af5
    cases.to_csv(path, sep="\t", index=False)
    r = subprocess.run(cmd, cwd=HERE, env=env, capture_output=True, text=True)
    if r.returncode == 0 or "only in v6" not in r.stdout:
        problems.append("check_af5.py did not report a changed frequency")
    path.write_text(original)
    return problems


def check(root: Path, diag_af: dict) -> list[str]:
    problems = check_bands() + check_identity(root)
    for study in ("snv", "replication"):
        refs = {"main": run(root, study, "main"), "bg1pct": run(root, study, "bg1pct")}
        res = run(root, study, "sweep")
        rows = res["rows"]
        if list(rows) != [f"{af:g}" for af in common.SWEEP_AF]:
            problems.append(f"{study}: thresholds {list(rows)}")
        for af, name in ((common.MAIN_AF, "main"), (common.LOOSE_AF, "bg1pct")):
            r = res["reproduction"][f"{af:g}"]
            if not r["identical"] or r["compared"] < 100:
                problems.append(f"{study} {af:g}: sweep's own check {r}")
            # The same check done here: every score and contrast both compute, in all cases and the gnomAD strata.
            mine, ref = leaves({k: rows[f"{af:g}"][k] for k in ("strata", "gnomad_split", "median_variants_per_case")}), \
                leaves({k: refs[name][k] for k in ("strata", "gnomad_split", "median_variants_per_case")})
            bad = [p for p, v in mine.items() if "verdict" not in p and (p not in ref or ref[p] != v)]
            if bad:
                problems.append(f"{study} {af:g}: {len(bad)} values differ from {name}, e.g. {bad[:3]}")
        medians = []
        for af in common.SWEEP_AF:
            row = rows[f"{af:g}"]
            lost = sum(x >= af for x in diag_af[study])
            if row["filter"]["cases_lost_diagnostic"] != lost or row["cases"] != len(diag_af[study]) - lost:
                problems.append(f"{study} {af:g}: {row['cases']} cases, {row['filter']['cases_lost_diagnostic']} "
                                f"lost; expected {lost} lost")
            has_evo2 = cb.EVO2 in row["strata"]["all"]["scores"]
            if has_evo2 != (af <= common.LOOSE_AF):
                problems.append(f"{study} {af:g}: Evo 2 {'reported' if has_evo2 else 'missing'}")
            for name in ("Rule + AVI minus Rule-based",
                         "Lookups + AVI (locked) minus Lookups combined (locked)" if study == "replication"
                         else "Lookups + AVI minus Lookups combined"):
                if name not in row["strata"]["all"]["contrasts"]:
                    problems.append(f"{study} {af:g}: no contrast {name}")
            if set(row["why"]["above_by_band"]) != {"AVI", "CADD"}:
                problems.append(f"{study} {af:g}: why table {row['why']['above_by_band']}")
            medians.append(row["median_variants_per_case"])
        if medians != sorted(medians):
            problems.append(f"{study}: median list shrinks as the threshold rises: {medians}")
        band = rows["0.05"]["why"]["background_by_band"]["1 to 5%"]
        if not band > 0 or rows["0.01"]["why"]["background_by_band"]["1 to 5%"] != 0:
            problems.append(f"{study}: 1 to 5% band {band} at 5%, should be empty at 1%")
        rho = rows["0.05"]["why"]["spearman_with_grpmax"]["AVI"]["rho"]
        if not rho < 0:
            problems.append(f"{study}: AVI should fall with grpmax in the synthetic data, rho {rho}")
    md = (root / "results" / "analysis_sweep.md").read_text()
    for af in common.SWEEP_AF:
        if f"| {af * 100:g}% |" not in md:
            problems.append(f"analysis_sweep.md has no row for {af * 100:g}%")
    if "—" in md or "–" in md:
        problems.append("analysis_sweep.md has an em or en dash")
    if "NOT identical" in md:
        problems.append("analysis_sweep.md reports a reproduction failure")
    env = os.environ | {"RDV2_RESULTS": str(root / "results")}
    subprocess.run([sys.executable, str(HERE / "figures" / "slide_charts.py"), "--brand", "gel", "--out",
                    str(root / "figures"), "--only", "af_sweep", "af_sweep_gain"], cwd=HERE, env=env, check=True,
                   stdout=subprocess.DEVNULL)
    for name in ("af_sweep_gel", "af_sweep_gain_gel"):
        for ext in ("png", "svg"):
            if not (root / "figures" / f"{name}.{ext}").stat().st_size:
                problems.append(f"figure {name}.{ext} is empty")
    return problems


def test_sweep() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        problems = check(root, synthetic(root))
    assert not problems, "\n".join(problems[:20])


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        problems = check(root, synthetic(root))
    for line in problems:
        print("  " + line)
    print(f"sweep: {len(problems)} problems")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
