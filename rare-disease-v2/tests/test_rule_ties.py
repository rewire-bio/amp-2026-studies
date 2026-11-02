"""The rule-based tie rule (PROTOCOL.md section 7) on hand-made cases, one case per branch.

Each case holds one diagnostic and one background variant; the test checks the background variant's outcome
(LOSES, TIED or BEATS) under Rule-based and Rule + AVI, and that competitors() and tied() count it.

  cd rare-disease-v2 && uv run -q --python 3.12 --with-requirements requirements.txt python tests/test_rule_ties.py

Exit status 0 when every case gives the expected outcome.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import case_benchmark as cb  # noqa: E402

N = np.nan
L, T, B = cb.LOSES, cb.TIED, cb.BEATS
COLS = ["var_class", "is_snv", "REVEL", "AlphaMissense", "SpliceAI", "CADD", "AVI"]

# name: (diagnostic, background, expected Rule-based, expected Rule + AVI)
CASES = {
    "more points": (("missense", True, 0.70, N, 0.0, 20, 1.0), ("nonsense", True, N, N, 0.0, 30, 0.0), B, B),
    "fewer points": (("nonsense", True, N, N, 0.0, 30, 0.0), ("missense", True, 0.95, N, 0.0, 30, 2.0), L, L),
    "same tool, higher score": (("missense", True, 0.70, N, 0.0, 30, 2.0), ("missense", True, 0.75, N, 0.0, 10, 1.0),
                                B, L),
    "same tool, lower score": (("missense", True, 0.75, N, 0.0, 10, 1.0), ("missense", True, 0.70, N, 0.0, 30, 2.0),
                               L, B),
    "same tool, equal score, CADD equal": (("intronic", True, N, N, 0.5, 20, 1.0),
                                           ("utr", True, N, N, 0.5, 20, 1.0), T, T),
    "same tool, equal score, CADD lower": (("intronic", True, N, N, 0.5, 20, 1.0),
                                           ("utr", True, N, N, 0.5, 19, 1.0), L, L),
    "different tools, CADD higher": (("intronic", True, N, N, 0.8, 20, N), ("missense", True, 0.70, N, 0.0, 25, N),
                                     B, B),
    "different tools, CADD lower": (("intronic", True, N, N, 0.8, 20, N), ("missense", True, 0.70, N, 0.0, 15, N),
                                    L, L),
    "REVEL against AlphaMissense: CADD decides": (("missense", True, 0.70, N, 0.0, 20, N),
                                                  ("missense", True, N, 0.80, 0.0, 20, N), T, T),
    "missense with SpliceAI strictly higher uses SpliceAI": (("intronic", True, N, N, 0.3, 30, N),
                                                             ("missense", True, 0.5, N, 0.4, 10, N), B, B),
    "no tool (PVS1), CADD equal": (("nonsense", True, N, N, N, 35, N), ("canonical_splice", True, N, N, 0.9, 35, N),
                                   T, T),
    "missing CADD counts against, not a tie": (("intronic", True, N, N, 0.8, 20, N),
                                               ("missense", True, 0.70, N, 0.0, N, N), B, B),
    "background with no points": (("intronic", True, N, N, 0.0, 5, N), ("utr", True, N, N, N, 40, N), L, L),
    "diagnostic with no points": (("intronic", True, N, N, N, 40, N), ("utr", True, N, N, 0.0, 5, N), B, B),
    "indels, SpliceAI higher": (("utr", False, N, N, 0.3, N, N), ("intronic", False, N, N, 0.5, N, N), B, B),
    "indels, SpliceAI equal": (("inframe", False, N, N, 0.0, N, N), ("inframe", False, N, N, 0.0, N, N), T, T),
    "indels, inframe without SpliceAI keeps +2": (("inframe", False, N, N, N, N, N),
                                                   ("inframe", False, N, N, N, N, N), T, T),
    "indels, scored beats missing": (("frameshift", False, N, N, N, N, N), ("frameshift", False, N, N, 0.0, N, N),
                                     B, B),
    "AVI scored beats AVI missing": (("intronic", True, N, N, 0.5, 30, N), ("utr", True, N, N, 0.5, 10, -5.0), L, B),
}


def frame() -> pd.DataFrame:
    rows = []
    for i, (name, (d, b, _, _)) in enumerate(CASES.items()):
        for diag, values in ((1, d), (0, b)):
            rows.append({"case": f"c{i:02d}", "diagnostic": diag, **dict(zip(COLS, values))})
    out = pd.DataFrame(rows)
    out["is_snv"] = out["is_snv"].astype(bool)
    return out


def test_rule_ties() -> list[str]:
    cases = frame()
    problems = []
    for n, col in ((cb.RULE, 2), (cb.RULE_AVI, 3)):
        out = cb.rule_outcomes(cases, with_avi=n == cb.RULE_AVI)
        cases[cb.outcome_column(n)] = out
        comp, tie = cb.competitors(cases, n), cb.tied(cases, n)
        for i, (name, spec) in enumerate(CASES.items()):
            got = out[2 * i + 1]
            if got != spec[col]:
                problems.append(f"{n}, {name}: expected {spec[col]}, got {got}")
            if out[2 * i] != cb.LOSES:
                problems.append(f"{n}, {name}: the diagnostic row is not LOSES")
            case = f"c{i:02d}"
            if comp[case] != int(spec[col] >= cb.TIED) or bool(tie[case]) != (spec[col] == cb.TIED):
                problems.append(f"{n}, {name}: competitors or tied do not match the outcome")
    assert not problems, "\n".join(problems)
    return problems


if __name__ == "__main__":
    try:
        test_rule_ties()
    except AssertionError as e:
        print(e)
        sys.exit(1)
    print(f"{len(CASES)} cases, Rule-based and Rule + AVI: all as expected")
