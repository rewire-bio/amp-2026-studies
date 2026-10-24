"""Compare the committed result files with an earlier tag: which point estimates and intervals changed.

Used when AVI was added (24 October 2026) to check that adding a score moved no point estimate, only the bootstrap
intervals (the draws are taken in score order, so a new score shifts the draws for the scores after it).

  python3 rare-disease/compare_results.py amp-2026-v1    # writes results/compare_<tag>.txt
"""
import json
import subprocess
import sys

from pathlib import Path

HERE = Path(__file__).resolve().parent
FILES = ["results.json", "results_indel.json", "results_cancer.json", "results_replication.json"]


def main() -> None:
    tag = sys.argv[1]
    lines = [f"Changes in rare-disease/results since {tag}"]
    for f in FILES:
        old = json.loads(subprocess.check_output(["git", "show", f"{tag}:rare-disease/results/{f}"], cwd=HERE))
        new = json.loads((HERE / "results" / f).read_text())
        points, intervals = [], 0
        for stratum, s in old["strata"].items():
            for section in ("scores", "added_to_lookups", "added_to_locked_lookups"):
                for name, v in s.get(section, {}).items():
                    nv = new["strata"][stratum][section][name]
                    for k in v:
                        if v[k] != nv[k]:
                            if k.endswith("_ci95"):
                                intervals += 1
                            else:
                                points.append(f"  {stratum} / {section} / {name} / {k}: {v[k]} -> {nv[k]}")
        added = sorted(set(new["coverage_diagnostic"]) - set(old["coverage_diagnostic"]))
        lines += [f"{f}: new scores {added or 'none'}; {intervals} intervals changed; "
                  f"{len(points)} point estimates changed"] + points
    out = HERE / "results" / f"compare_{tag}.txt"
    out.write_text("\n".join(lines) + "\n")
    print(out.read_text())


if __name__ == "__main__":
    main()
