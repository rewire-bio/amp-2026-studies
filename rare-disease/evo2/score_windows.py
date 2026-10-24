"""Score the rare-disease windows with Evo 2 7B on the pod. Default numerics, batch size 1, one GPU type.

Writes ref_scores.tsv and var_scores.tsv as it goes (resumable), then scores.tsv with
delta = variant mean log-likelihood minus reference mean log-likelihood, and receipt.json.

  python score_windows.py --data /work/data --out /work/out
"""
import argparse
import csv
import datetime as dt
import gzip
import json
import sys
import time
import traceback
from pathlib import Path

import torch

sys.path.insert(0, "/work/code/src")
sys.path.insert(0, "/work/code/followup")
from gpu_followup import load_model  # noqa: E402

CHUNK = 64


def done_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open() as h:
        return {line.split("\t", 1)[0] for line in h if line.strip()}


def score_file(model, rows: list[tuple[str, str]], out: Path, label: str) -> tuple[int, float]:
    """Score (id, seq) rows not already in out; append id and score per line."""
    seen = done_ids(out)
    todo = [r for r in rows if r[0] not in seen]
    start = time.perf_counter()
    with out.open("a") as h:
        for i in range(0, len(todo), CHUNK):
            chunk = todo[i:i + CHUNK]
            scores = model.score_sequences([s for _, s in chunk], batch_size=1, reduce_method="mean")
            for (key, _), s in zip(chunk, scores):
                h.write(f"{key}\t{float(s)!r}\n")
            h.flush()
            print(f"{label}: {len(seen) + i + len(chunk)}/{len(rows)}", flush=True)
    torch.cuda.synchronize()
    return len(todo), time.perf_counter() - start


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("/work/data"))
    p.add_argument("--out", type=Path, default=Path("/work/out"))
    p.add_argument("--cache", type=Path, default=Path("/work/hf"))
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    receipt = {"started_utc": dt.datetime.now(dt.UTC).isoformat(), "gpu_name": torch.cuda.get_device_name(0),
               "status": "running"}
    try:
        t = time.perf_counter()
        model, applied = load_model("default", a.cache)
        receipt.update(applied=applied, load_s=round(time.perf_counter() - t, 3))
        with gzip.open(a.data / "refs.tsv.gz", "rt") as h:
            refs = [(r["ref_id"], r["seq"]) for r in csv.DictReader(h, delimiter="\t")]
        with gzip.open(a.data / "vars.tsv.gz", "rt") as h:
            variants = list(csv.DictReader(h, delimiter="\t"))
        var_rows = [(f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}", v["seq"]) for v in variants]
        n_ref, ref_s = score_file(model, refs, a.out / "ref_scores.tsv", "ref")
        n_var, var_s = score_file(model, var_rows, a.out / "var_scores.tsv", "var")

        def load(path):
            with path.open() as h:
                return {k: float(v) for k, v in (line.rstrip("\n").split("\t") for line in h if line.strip())}

        ref_scores, var_scores = load(a.out / "ref_scores.tsv"), load(a.out / "var_scores.tsv")
        with (a.out / "scores.tsv").open("w", newline="") as h:
            w = csv.writer(h, delimiter="\t")
            w.writerow(["chrom", "pos", "ref", "alt", "ref_score", "var_score", "evo2_delta"])
            for v, (key, _) in zip(variants, var_rows):
                r, s = ref_scores[v["ref_id"]], var_scores[key]
                w.writerow([v["chrom"], v["pos"], v["ref"], v["alt"], repr(r), repr(s), repr(s - r)])
        windows = n_ref + n_var
        seconds = ref_s + var_s
        receipt.update(reference_windows=len(refs), variant_windows=len(var_rows), windows_scored_this_run=windows,
                       scoring_s=round(seconds, 1), windows_per_s=round(windows / seconds, 4) if seconds else None,
                       peak_memory_gib=round(torch.cuda.max_memory_allocated() / 1024 ** 3, 2), status="completed")
    except Exception:
        receipt.update(status="failed", error=traceback.format_exc())
    finally:
        receipt["finished_utc"] = dt.datetime.now(dt.UTC).isoformat()
        (a.out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
