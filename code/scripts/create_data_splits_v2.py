"""GRIMM v2 splitting — paper-faithful, leakage-free.

Fixes vs v1 (create_data_splits_custom.py):
  Bug A: test-2 is the held-out orphan set only (never also written to train).
  Bug B: low-support labels (1-2 clusters) are split BY CLUSTER, not by row.
  Bug D: cluster->partition assignment is shuffled per fold (folds are independent).
  Bug F: all randomness is seeded (reproducible).
(Bug C — real per-protein sequences — is handled when sequences are joined in a
later step; splitting itself only needs accession + UniRef50.)

Per-label rules (cluster counts), matching the preprint:
  >=10 clusters : ceil(10%) test1, ceil(10%) valid, rest train   (~80/10/10)
  6-9  clusters : 2 test1, 2 valid, rest train
  3-5  clusters : 1 test1, 1 valid, rest train
  2    clusters : 1 test1, 0 valid, 1 train          (by cluster — Bug B fix)
  1    cluster  : orphan -> ~20% test2 / ~80% train across the 5 folds (no reuse)

Cluster-holdout is enforced PER LABEL: a cluster shared by two labels may sit in
different partitions for each label (per-label stratification, as in the paper),
but within a label no cluster is ever split across partitions.

Output: data/v2/split_{i}/assignments.tsv  (Entry, EC number, UniRef50, partition)
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

LABEL = "EC number"
CLUSTER = "UniRef50"
SEED = 10297362
N_SPLITS = 5


def partition_sizes(n: int) -> tuple[int, int, int]:
    """(n_test1, n_valid, n_train) cluster counts for a label with n>=2 clusters."""
    if n >= 10:
        n_test = math.ceil(0.1 * n)
        n_valid = math.ceil(0.1 * n)
    elif n >= 6:
        n_test, n_valid = 2, 2
    elif n >= 3:
        n_test, n_valid = 1, 1
    else:  # n == 2
        n_test, n_valid = 1, 0
    return n_test, n_valid, n - n_test - n_valid


def split(df: pd.DataFrame, out_root: Path) -> None:
    # cluster -> partition is decided per (label); collect orphan clusters globally
    # orphans are deterministically chunked across folds (seeded), so each orphan
    # cluster is test2 in exactly one fold and train in the others.
    rng = np.random.RandomState(SEED)
    label_clusters = df.groupby(LABEL)[CLUSTER].unique()

    orphan_clusters = []  # (label, cluster) with exactly 1 cluster for that label
    splittable = {}       # label -> np.array of its clusters (>=2)
    for label, clusters in label_clusters.items():
        if len(clusters) == 1:
            orphan_clusters.append((label, clusters[0]))
        else:
            splittable[label] = clusters

    orphan_clusters = np.array(orphan_clusters, dtype=object)
    perm = rng.permutation(len(orphan_clusters))
    orphan_clusters = orphan_clusters[perm]
    orphan_folds = np.array_split(orphan_clusters, N_SPLITS)  # chunk i -> test2 in fold i

    for i in range(N_SPLITS):
        fold_rng = np.random.RandomState(SEED + i)
        # (label, cluster) -> partition
        assign: dict[tuple, str] = {}

        for label, clusters in splittable.items():
            clusters = np.array(clusters)
            fold_rng.shuffle(clusters)  # Bug D: shuffle which cluster goes where
            n_test, n_valid, n_train = partition_sizes(len(clusters))
            for j, c in enumerate(clusters):
                if j < n_test:
                    part = "test1"
                elif j < n_test + n_valid:
                    part = "valid"
                else:
                    part = "train"
                assign[(label, c)] = part

        test2_set = {tuple(x) for x in orphan_folds[i]}
        for k in range(N_SPLITS):
            for lc in orphan_folds[k]:
                key = tuple(lc)
                assign[key] = "test2" if k == i else "train"

        df_fold = df.copy()
        df_fold["partition"] = [
            assign[(row_label, row_cluster)]
            for row_label, row_cluster in zip(df_fold[LABEL], df_fold[CLUSTER])
        ]

        out_dir = out_root / f"split_{i + 1}"
        out_dir.mkdir(parents=True, exist_ok=True)
        cols = ["Entry", LABEL, CLUSTER, "partition"]
        df_fold[cols].to_csv(out_dir / "assignments.tsv", sep="\t", index=False)
        counts = df_fold["partition"].value_counts()
        print(f"split {i+1}: " + "  ".join(
            f"{p}={int(counts.get(p, 0))}" for p in ["train", "valid", "test1", "test2"]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base-table", type=Path, default=Path("v2_build/data/base_table.tsv"))
    p.add_argument("--out-root", type=Path, default=Path("v2_build/data/v2_splits"))
    args = p.parse_args()
    df = pd.read_csv(args.base_table, sep="\t")
    print(f"base table: {len(df)} accessions, {df[LABEL].nunique()} labels, "
          f"{df[CLUSTER].nunique()} clusters")
    split(df, args.out_root)


if __name__ == "__main__":
    main()
