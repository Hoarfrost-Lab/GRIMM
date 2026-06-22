"""Decontaminate the deployed GRIMM v1 test-2 splits.

Bug A in create_data_splits_custom.py wrote the full `orphans` set to test2.csv
instead of the held-out `test2` sample, so ~80% of test-2 rows are also in train.
This script repairs the DEPLOYED files directly (we do NOT regenerate v1 — see
Bug F / unseeded RNG): for each fold it drops any test-2 row whose accession OR
sequence appears in that fold's train set, leaving a true out-of-distribution set.

Exact-match definition mirrors teacup/evaluation/decontaminate.py (reviewer-defensible).

Usage:
    python code/scripts/decontaminate_test2.py \
        --in-root  /path/to/hf_dataset/EC \
        --out-root /path/to/patched/EC
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# train/test file names differ slightly by modality (valid vs validation)
MODALITIES = {
    "amino_acids": {"train": "train.csv", "test2": "test2.csv",
                    "test1": "test1.csv", "valid": "validation.csv"},
    "nucleotides": {"train": "train.csv", "test2": "test2.csv",
                    "test1": "test1.csv", "valid": "valid.csv"},
}


def decontaminate_fold(split_dir: Path, names: dict) -> tuple[pd.DataFrame, dict]:
    train = pd.read_csv(split_dir / names["train"], sep="\t")
    valid = pd.read_csv(split_dir / names["valid"], sep="\t")
    test1 = pd.read_csv(split_dir / names["test1"], sep="\t")
    test2 = pd.read_csv(split_dir / names["test2"], sep="\t")

    # Fully disjoint test-2: exclude any row whose accession OR sequence appears in
    # ANY other partition (train, valid, test-1). Reviewer-defensible exact-match OOD.
    other_entries = set(train["Entry"]) | set(valid["Entry"]) | set(test1["Entry"])
    other_seqs = set(train["Sequence"]) | set(valid["Sequence"]) | set(test1["Sequence"])

    keep = ~(test2["Entry"].isin(other_entries) | test2["Sequence"].isin(other_seqs))
    patched = test2[keep].copy()

    stats = {
        "n_before": len(test2),
        "n_after": len(patched),
        "dropped": len(test2) - len(patched),
        "resid_entry": int(patched["Entry"].isin(other_entries).sum()),
        "resid_seq": int(patched["Sequence"].isin(other_seqs).sum()),
    }
    return patched, stats


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in-root", type=Path, required=True,
                   help="EC dir containing amino_acids/ and nucleotides/")
    p.add_argument("--out-root", type=Path, required=True)
    p.add_argument("--splits", type=int, default=5)
    args = p.parse_args()

    print(f"{'modality':12} {'split':7} {'before':>7} {'after':>7} {'dropped':>8} "
          f"{'resid_E':>8} {'resid_S':>8}")
    print("-" * 56)
    for modality, names in MODALITIES.items():
        for s in range(1, args.splits + 1):
            split_dir = args.in_root / modality / f"split_{s}"
            if not split_dir.exists():
                print(f"{modality:12} split_{s}: MISSING ({split_dir})")
                continue
            patched, st = decontaminate_fold(split_dir, names)
            out_dir = args.out_root / modality / f"split_{s}"
            out_dir.mkdir(parents=True, exist_ok=True)
            patched.to_csv(out_dir / "test2.csv", sep="\t", index=False)
            flag = "" if (st["resid_entry"] == 0 and st["resid_seq"] == 0) else "  <-- LEAK!"
            print(f"{modality:12} split_{s} {st['n_before']:7d} {st['n_after']:7d} "
                  f"{st['dropped']:8d} {st['resid_entry']:8d} {st['resid_seq']:8d}{flag}")


if __name__ == "__main__":
    main()
