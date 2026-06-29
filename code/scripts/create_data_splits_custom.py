"""GRIMM-EC amino-acid splits (v2) — paper-faithful, leakage-free.

Per-label (EC) stratification by UniRef50 cluster. A cluster shared by two labels
may sit in different partitions for each label (per-label stratification, as in the
preprint), but within a label no cluster is ever split across partitions.

Rules by cluster count per label (matching the preprint):
  >=10 clusters : ceil(10%) test1, ceil(10%) valid, rest train   (~80/10/10)
  6-9  clusters : 2 test1, 2 valid, rest train
  3-5  clusters : 1 test1, 1 valid, rest train
  2    clusters : 1 test1, 0 valid, 1 train      (by cluster)
  1    cluster  : orphan -> ~20% test2 / ~80% train across the 5 folds (no reuse)

test2 is the held-out open-set (labels absent from training). All randomness is
seeded; cluster->partition assignment is shuffled per fold (folds are independent).
A final (Sequence, EC) guard drops any eval row whose exact sequence+label is in
train (preserves intended cross-label homology).

This reproduces the method described in the README/preprint. (The pre-v2 history of
this file in git contains the original implementation, which had the leakage bugs.)
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

LABEL, CLUSTER, SEED, N_SPLITS = "EC number", "UniRef50", 10297362, 5
COLS = ["Entry", "EC number", "EMBL", "RefSeq", "UniRef50", "UniRef90", "UniRef100", "Sequence"]
PART_FILE = {"train": "train.tsv", "valid": "validation.tsv",
             "test1": "test1.tsv", "test2": "test2.tsv"}


def partition_sizes(n: int) -> tuple:
    if n >= 10:
        return math.ceil(0.1 * n), math.ceil(0.1 * n), None
    if n >= 6:
        return 2, 2, None
    if n >= 3:
        return 1, 1, None
    return 1, 0, None  # n == 2


def assign_fold(df, splittable, orphan_folds, i):
    fold_rng = np.random.RandomState(SEED + i)
    assign = {}
    for label, clusters in splittable.items():
        clusters = np.array(clusters)
        fold_rng.shuffle(clusters)
        n_test, n_valid, _ = partition_sizes(len(clusters))
        for j, c in enumerate(clusters):
            assign[(label, c)] = ("test1" if j < n_test else
                                  "valid" if j < n_test + n_valid else "train")
    for k in range(N_SPLITS):
        for label, c in orphan_folds[k]:
            assign[(label, c)] = "test2" if k == i else "train"
    return [assign[(lab, cl)] for lab, cl in zip(df[LABEL], df[CLUSTER])]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base-table", type=Path, default=Path("v2_build/data/base_table_with_seq.tsv"))
    p.add_argument("--out-root", type=Path, default=Path("v2_build/data/GRIMM_v2/amino_acids"))
    args = p.parse_args()

    df = pd.read_csv(args.base_table, sep="\t")
    n_missing = df["Sequence"].isna().sum()
    df = df[df["Sequence"].notna()].reset_index(drop=True)
    print(f"{len(df)} accessions with sequence ({n_missing} dropped); "
          f"{df[LABEL].nunique()} labels, {df[CLUSTER].nunique()} clusters")

    label_clusters = df.groupby(LABEL)[CLUSTER].unique()
    orphans = np.array([(lab, cs[0]) for lab, cs in label_clusters.items() if len(cs) == 1], dtype=object)
    splittable = {lab: cs for lab, cs in label_clusters.items() if len(cs) > 1}
    orphans = orphans[np.random.RandomState(SEED).permutation(len(orphans))]
    orphan_folds = np.array_split(orphans, N_SPLITS)

    print(f"\n{'split':7} {'train':>8} {'valid':>8} {'test1':>8} {'test2':>8}  {'t1_leak':>8} {'t2_leak':>8}")
    for i in range(N_SPLITS):
        df["partition"] = assign_fold(df, splittable, orphan_folds, i)
        frames = {p: df[df["partition"] == p][COLS].copy() for p in PART_FILE}

        # (Sequence, EC) guard — drop eval rows whose exact sequence+label is in train.
        train_pairs = set(zip(frames["train"]["Sequence"], frames["train"]["EC number"]))
        for part in ("valid", "test1", "test2"):
            f = frames[part]
            frames[part] = f[[(s, e) not in train_pairs for s, e in zip(f["Sequence"], f["EC number"])]]

        out_dir = args.out_root / f"split_{i+1}"
        out_dir.mkdir(parents=True, exist_ok=True)
        for part, fname in PART_FILE.items():
            frames[part].to_csv(out_dir / fname, sep="\t", index=False)

        t1 = sum((s, e) in train_pairs for s, e in zip(frames["test1"]["Sequence"], frames["test1"]["EC number"]))
        t2 = sum((s, e) in train_pairs for s, e in zip(frames["test2"]["Sequence"], frames["test2"]["EC number"]))
        print(f"split_{i+1} {len(frames['train']):8d} {len(frames['valid']):8d} "
              f"{len(frames['test1']):8d} {len(frames['test2']):8d}  {t1:8d} {t2:8d}")


if __name__ == "__main__":
    main()
