#!/usr/bin/env python3
"""Summarise record_run.sh's scratch runs into run/: summary.json, ranks.tsv and trace_first.txt.

  python3 summarise_run.py --work <record_run.sh scratch dir> --ops <evo2-cloud-operations> --out run

Called by record_run.sh. Standard library only. CADD and SpliceAI values are read to compute ranks and ties but are
not written out. Ranks are 1 for the highest score; a tied score has a range of ranks (min to max), and the pipeline's
own rank inside a tie follows the VCF order (Python's sort is stable).
"""
import argparse
import csv
import gzip
import json
import re
import shutil
import statistics
from collections import Counter
from pathlib import Path

SCORES = {"cadd": "cadd_raw", "spliceai": "spliceai_max", "evo2": "evo2_lof"}


def ranked(path: Path) -> dict[str, dict]:
    with path.open() as handle:
        return {r["variant"]: r for r in csv.DictReader(handle, delimiter="\t")}


def tie_range(values: dict[str, float], key: str) -> tuple[int, int]:
    v = values[key]
    return 1 + sum(x > v for x in values.values()), sum(x >= v for x in values.values())


def vcf_body(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [l for l in handle if not l.startswith("##bcftools")]


def spearman(a: list[float], b: list[float]) -> float:
    def ranks(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    ra, rb = ranks(a), ranks(b)
    ma, mb = statistics.fmean(ra), statistics.fmean(rb)
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    return cov / (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5


def text(path: Path) -> str:
    return path.read_text().strip() if path.exists() else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--ops", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    w, ops = args.work, args.ops
    args.out.mkdir(parents=True, exist_ok=True)
    pipe = Path(__file__).resolve().parents[1]

    truth = (pipe / "test/case_brca1.truth.txt").read_text().split("\t")[0]
    first = ranked(w / "out_first/ranked.tsv")
    ids = list(first)
    values = {s: {k: float(r[c]) for k, r in first.items() if r[c] != "missing"} for s, c in SCORES.items()}
    l40s = ranked(w / "out_l40s/ranked.tsv")
    p3 = ranked(w / "out_p3/ranked.tsv")
    def full(name: str) -> dict[str, float]:  # unrounded Evo 2 scores from the step output (ranked.tsv has 6 decimals)
        with (w / f"out_{name}/evo2_scores.tsv").open() as handle:
            return {f"{r['chrom']}:{r['pos']}:{r['ref']}>{r['alt']}": float(r["evo2_lof"])
                    for r in csv.DictReader(handle, delimiter="\t")}
    l40s_values, h100_values = full("l40s"), full("first")

    # Case composition from the study's variant table (aggregate counts only)
    with (ops / "data/derived/brca1_variants.tsv").open() as handle:
        table = {f"17:{v['pos']}:{v['ref']}>{v['alt']}": v for v in csv.DictReader(handle, delimiter="\t")}
    case = [table[k] for k in ids]
    positions = len({v["pos"] for v in case})

    # Replayed scores must equal the recorded H100 scores (reference minus variant, the recorded delta negated)
    with (ops / "evidence/runs/recorded-b75c67aa/platforms/P2/attempt-1/w1_brca1_scores.tsv").open() as handle:
        recorded = {r["id"].removeprefix("chr"): -float(r["delta"]) for r in csv.DictReader(handle, delimiter="\t")}
    replay_max_diff = max(abs(h100_values[k] - recorded[k]) for k in ids)
    values["evo2"] = h100_values

    # Write ranks.tsv
    with (args.out / "ranks.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant", "diagnostic", "consequence", "cadd_rank", "cadd_rank_min", "cadd_rank_max",
                         "spliceai_rank", "spliceai_rank_min", "spliceai_rank_max", "evo2_lof_h100", "evo2_rank_h100",
                         "evo2_lof_l40s", "evo2_rank_l40s"])
        for k in ids:
            r = first[k]
            writer.writerow([k, int(k == truth), table[k]["consequence"],
                             r["cadd_rank"], *tie_range(values["cadd"], k),
                             r["spliceai_rank"], *tie_range(values["spliceai"], k),
                             r["evo2_lof"], r["evo2_rank"], l40s[k]["evo2_lof"], l40s[k]["evo2_rank"]])

    def diag(scores: dict[str, float], rank: str) -> dict:
        lo, hi = tie_range(scores, truth)
        return {"rank": int(rank), "rank_min": lo, "rank_max": hi, "variants_sharing_the_score": hi - lo + 1}

    diagnostic = {
        "variant": truth, "consequence": table[truth]["consequence"],
        "cadd": diag(values["cadd"], first[truth]["cadd_rank"]),
        "spliceai": diag(values["spliceai"], first[truth]["spliceai_rank"]),
        "evo2_h100": diag(values["evo2"], first[truth]["evo2_rank"]),
        "evo2_l40s": diag(l40s_values, l40s[truth]["evo2_rank"]),
        "evo2_lof_h100": h100_values[truth], "evo2_lof_l40s": l40s_values[truth],
    }

    # Determinism, caching and the two guards
    second_same_ranked = (w / "out_first/ranked.tsv").read_bytes() == (w / "out_second/ranked.tsv").read_bytes()
    annotated_bytes_same = (w / "out_first/annotated.vcf.gz").read_bytes() == (w / "out_second/annotated.vcf.gz").read_bytes()
    annotated_body_same = vcf_body(w / "out_first/annotated.vcf.gz") == vcf_body(w / "out_second/annotated.vcf.gz")
    cached = len(re.findall(r"Cached process", text(w / "log_resume.txt")))
    submitted_resume = len(re.findall(r"Submitted process", text(w / "log_resume.txt")))
    diffs = [abs(h100_values[k] - l40s_values[k]) for k in ids]
    rank_moves = [abs(int(first[k]["evo2_rank"]) - int(l40s[k]["evo2_rank"])) for k in ids]
    biggest = ids[rank_moves.index(max(rank_moves))]
    top = lambda table, n: {k for k in ids if int(table[k]["evo2_rank"]) <= n}  # noqa: E731
    norm_stats = re.search(r"total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\s*([\d/]+)",
                           text(w / "err_mismatch.txt"))
    mismatch_rows = len(ranked(w / "out_mismatch/ranked.tsv")) if (w / "out_mismatch/ranked.tsv").exists() else None
    missing_msg = re.search(r"No recorded Evo 2 score for \S+", text(w / "log_missing.txt") + text(w / "err_missing.txt"))

    trace = list(csv.DictReader((w / "trace_first.txt").open(), delimiter="\t"))
    shutil.copy(w / "trace_first.txt", args.out / "trace_first.txt")

    ops_run = json.loads((ops / "evidence/runs/recorded-b75c67aa/operations.json").read_text())["platforms"]
    agreement = json.loads((ops / "evidence/runs/recorded-b75c67aa/agreement-details.json").read_text())
    cloud = {}
    for pid in ("P2", "P3", "P5"):
        p = ops_run[pid]
        cloud[pid] = {"gpu": p["gpu"], "provider": p["instance"]["provider"],
                      "purchase_option": p["instance"]["purchase_option"],
                      "instance_type": p["instance"]["instance_type"], "region": p["instance"]["region"],
                      "attempt": p["attempt"], "cost": p["cost"], "docker_pull_s": p["docker_pull_s"],
                      "events_utc": p["events_utc"], "timings_s": p["timings_s"],
                      "w1_sequences_per_second": p["w1_sequences_per_second"],
                      "w1_peak_memory_gib": p["w1_peak_memory_gib"]}
    rate = ops_run["P2"]["w1_sequences_per_second"]
    sequences = positions + len(ids)

    summary = {
        "versions": dict(line.split("\t", 1) for line in text(w / "versions.tsv").splitlines()),
        "test_case_rebuilt_identical": text(w / "test_check.txt") == "identical",
        "windows_match_source_study": json.loads(text(w / "windows_check.json")),
        "case": {"variants": len(ids), "positions": positions, "evo2_sequences": sequences,
                 "consequences": dict(Counter(v["consequence"] for v in case).most_common()),
                 "assay_classes": dict(Counter(v["func_class"] for v in case).most_common())},
        "coverage": {s: len(values[s]) for s in SCORES},
        "diagnostic": diagnostic,
        "replay_matches_recorded_max_abs_diff": replay_max_diff,
        "runs": {name: {"exit": int(text(w / f"exit_{name}.txt")), "wall_s": float(text(w / f"wall_{name}.txt"))}
                 for name in ("first", "second", "resume", "p3", "l40s", "mismatch", "missing")},
        "first_run_tasks": [{"name": t["name"], "status": t["status"], "duration": t["duration"]} for t in trace],
        "determinism": {"ranked_identical": second_same_ranked, "annotated_vcf_bytes_identical": annotated_bytes_same,
                        "annotated_vcf_identical_without_bcftools_command_lines": annotated_body_same},
        "resume": {"cached_tasks": cached, "submitted_tasks": submitted_resume},
        "p3_ranked_identical": (w / "out_first/ranked.tsv").read_bytes() == (w / "out_p3/ranked.tsv").read_bytes(),
        "p3_evo2_identical": all(first[k]["evo2_lof"] == p3[k]["evo2_lof"] for k in ids),
        "l40s_vs_h100": {"n": len(ids), "median_abs_diff": statistics.median(diffs), "max_abs_diff": max(diffs),
                         "spearman": spearman([h100_values[k] for k in ids], [l40s_values[k] for k in ids]),
                         "variants_whose_evo2_rank_changes": sum(m > 0 for m in rank_moves),
                         "median_rank_move": statistics.median(rank_moves),
                         "variants_moving_at_most_2": sum(m <= 2 for m in rank_moves),
                         "variants_moving_5_or_more": sum(m >= 5 for m in rank_moves),
                         "variants_moving_10_or_more": sum(m >= 10 for m in rank_moves),
                         "largest_rank_move": max(rank_moves),
                         "largest_move_variant": biggest, "largest_move_h100_rank": int(first[biggest]["evo2_rank"]),
                         "largest_move_l40s_rank": int(l40s[biggest]["evo2_rank"]),
                         "top5_ranks_unchanged": all(first[k]["evo2_rank"] == l40s[k]["evo2_rank"]
                                                     for k in ids if int(first[k]["evo2_rank"]) <= 5),
                         "top10_shared_variants": len(top(first, 10) & top(l40s, 10))},
        "guards": {"ref_mismatch": {"exit": int(text(w / "exit_mismatch.txt")),
                                    "bcftools_norm_counts": norm_stats.group(1) if norm_stats else None,
                                    "ranked_rows": mismatch_rows},
                   "missing_score": {"exit": int(text(w / "exit_missing.txt")),
                                     "message": missing_msg.group(0) if missing_msg else None}},
        "cloud_run_recorded_b75c67aa": cloud,
        "cloud_agreement": agreement,
        "case_gpu_estimate": {"sequences": sequences, "rate_sequences_per_s_P2": rate,
                              "scoring_seconds": sequences / rate},
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("diagnostic", "determinism", "resume", "l40s_vs_h100", "guards")}, indent=1))


if __name__ == "__main__":
    main()
