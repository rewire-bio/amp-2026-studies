"""Adding or removing one score must not change any result that does not involve it (PROTOCOL.md sections 11, 12).

Runs case_benchmark.py on one study and analysis twice, with and without --drop-score (GPN-Star-M by default), and
checks that every result row that does not name the dropped score (in its path or its value) is identical in both
files, and that the run without it has no row the full run lacks. Rows are the leaves of the JSON (one value each,
intervals whole); in lists of score names the dropped name is removed before comparing. Only dropped_scores, which
names the dropped score by design, is skipped.

Uses common.DATA and common.RESULTS, so point RDV2_DATA and RDV2_RESULTS at a scratch copy to test there.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python tests/test_invariance.py \\
      [snv|indel|replication] [--analysis main] [--drop GPN-Star-M]

For added scores only (AVI, GPN-MSA, GPN-Star-V, GPN-Star-M, Evo 2). Dropping a lookup (CADD, SpliceAI, REVEL,
AlphaMissense) changes Lookups combined and the rule-based scores, whose names do not show it, so those rows differ
by design.

Exit status 0 when every row matches.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import common  # noqa: E402

SKIP = {"dropped_scores"}
ADDED = ["AVI", "GPN-MSA", "GPN-Star-V", "GPN-Star-M", "Evo 2"]


def rows(x, path=()) -> dict[tuple, object]:
    """Leaves of a JSON value, keyed by their path; lists are leaves (intervals, merge lists)."""
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            out |= rows(v, path + (k,))
        return out
    return {path: x}


def involves(path: tuple, value, name: str) -> bool:
    """A row involves the score if its path or its string value names it."""
    return any(name in str(p) for p in path) or (isinstance(value, str) and name in value)


def without_name(rows_: dict[tuple, object], name: str) -> dict[tuple, object]:
    """Drop the rows that involve the score; in lists of names (such as scores_without_coverage), drop the name."""
    out = {}
    for p, v in rows_.items():
        if p[0] in SKIP or involves(p, v, name):
            continue
        if isinstance(v, list):
            v = [x for x in v if not (isinstance(x, str) and name in x)]
        out[p] = v
    return out


def run(study: str, analysis: str, drop: str | None) -> dict:
    cmd = [sys.executable, str(HERE / "case_benchmark.py"), "--study", study, "--analysis", analysis]
    if drop:
        cmd += ["--drop-score", drop]
    subprocess.run(cmd, cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    suffix = f"_without_{re.sub('[^A-Za-z0-9]+', '_', drop)}" if drop else ""  # as case_benchmark.py names it
    return json.loads((common.RESULTS / f"{study}_{analysis}{suffix}.json").read_text())


def compare(study: str, analysis: str, drop: str) -> list[str]:
    full, without = rows(run(study, analysis, None)), rows(run(study, analysis, drop))
    kept, rest = without_name(full, drop), without_name(without, drop)
    problems = []
    for p, v in kept.items():
        if p not in rest:
            problems.append(f"missing without {drop}: {' / '.join(map(str, p))}")
        elif json.dumps(rest[p]) != json.dumps(v):
            problems.append(f"differs: {' / '.join(map(str, p))}: {v} vs {rest[p]}")
    for p, v in without.items():
        if p[0] not in SKIP and involves(p, v, drop):
            problems.append(f"still present without {drop}: {' / '.join(map(str, p))}")
    problems += [f"only in the run without {drop}: {' / '.join(map(str, p))}" for p in rest if p not in kept]
    print(f"{study} {analysis}: {len(kept)} rows compared, {len(full) - len(kept)} rows involve {drop} or are "
          f"skipped, {len(problems)} problems")
    return problems


def test_invariance(study: str = "snv", analysis: str = "main", drop: str = "GPN-Star-M") -> None:
    problems = compare(study, analysis, drop)
    assert not problems, "\n".join(problems[:20])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("study", nargs="?", default="snv", choices=["snv", "indel", "replication"])
    p.add_argument("--analysis", default="main")
    p.add_argument("--drop", default="GPN-Star-M")
    a = p.parse_args()
    if a.drop not in ADDED:
        p.error(f"--drop must be an added score: {', '.join(ADDED)}")
    problems = compare(a.study, a.analysis, a.drop)
    for line in problems[:40]:
        print("  " + line)
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
