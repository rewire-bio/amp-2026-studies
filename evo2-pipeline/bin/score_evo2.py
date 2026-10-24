#!/usr/bin/env python3
"""Score windows with Evo 2 7B on a GPU. Runs inside the pinned evo2-cloud-operations image.

The checkpoint is downloaded at a pinned revision and its SHA-256 checked before it is loaded, because
torch.load unpickles the file. Each reference window is scored once and shared by its variants.
"""
import argparse
import csv
import hashlib

WEIGHTS = {"repo_id": "arcinstitute/evo2_7b", "filename": "evo2_7b.pt",
           "revision": "bda0089f92582d5baabf0f22d9fc85f3588f6b58",
           "sha256": "c66645929dc1b9c631f5be656da8726f38946315dc9167000a615dd626fcecf4"}


def verified_checkpoint() -> str:
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(repo_id=WEIGHTS["repo_id"], filename=WEIGHTS["filename"], revision=WEIGHTS["revision"])
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    if digest.hexdigest() != WEIGHTS["sha256"]:
        raise SystemExit(f"SHA-256 mismatch for {path}")
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    from evo2 import Evo2

    model = Evo2("evo2_7b", local_path=verified_checkpoint())
    with open(args.windows) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ref_scores: dict[str, float] = {}
    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", "ref_ll", "var_ll", "evo2_lof"])
        for row in rows:
            if row["ref_window"] not in ref_scores:
                ref_scores[row["ref_window"]] = float(model.score_sequences([row["ref_seq"]], batch_size=1, reduce_method="mean")[0])
            ref_ll = ref_scores[row["ref_window"]]
            var_ll = float(model.score_sequences([row["var_seq"]], batch_size=1, reduce_method="mean")[0])
            writer.writerow([row["chrom"], row["pos"], row["ref"], row["alt"], repr(ref_ll), repr(var_ll), repr(ref_ll - var_ll)])


if __name__ == "__main__":
    main()
