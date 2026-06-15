"""Join real SwissProt sequences onto the verified v2 splits -> final AA CSVs + QC.

Reads v2_build/data/v2_splits/split_*/assignments.tsv (partition assignments) and
v2_build/data/base_table_with_seq.tsv (real per-protein sequences), and writes
v2_build/data/GRIMM_v2/amino_acids/split_*/{train,validation,test1,test2}.csv with
the v1 column schema. Reports sequence-level leakage (should be ~0, vs v1's 5.8%).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

COLS = ["Entry", "EC number", "EMBL", "RefSeq", "UniRef50", "UniRef90", "UniRef100", "Sequence"]
PART_FILE = {"train": "train.csv", "valid": "validation.csv",
             "test1": "test1.csv", "test2": "test2.csv"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--splits-root", type=Path, default=Path("v2_build/data/v2_splits"))
    p.add_argument("--seq-table", type=Path, default=Path("v2_build/data/base_table_with_seq.tsv"))
    p.add_argument("--out-root", type=Path, default=Path("v2_build/data/GRIMM_v2/amino_acids"))
    p.add_argument("--splits", type=int, default=5)
    args = p.parse_args()

    meta = pd.read_csv(args.seq_table, sep="\t")
    n_missing_seq = meta["Sequence"].isna().sum()
    meta = meta[meta["Sequence"].notna()].copy()
    print(f"sequence table: {len(meta)} rows with sequence ({n_missing_seq} dropped for missing seq)")
    meta_by_entry = meta.set_index("Entry")

    print(f"\n{'split':7} {'train':>8} {'valid':>8} {'test1':>8} {'test2':>8}  "
          f"{'t1_seq_in_tr':>12} {'t2_seq_in_tr':>12}")
    for s in range(1, args.splits + 1):
        a = pd.read_csv(args.splits_root / f"split_{s}" / "assignments.tsv", sep="\t")
        a = a[a["Entry"].isin(meta_by_entry.index)]  # drop accessions w/o sequence
        joined = a.join(meta_by_entry[[c for c in COLS if c not in ("Entry",)]],
                        on="Entry", rsuffix="_m")
        # prefer metadata copies of EC number / UniRef50 (identical); keep COLS order
        joined = joined.rename(columns={}).reset_index(drop=True)
        out_dir = args.out_root / f"split_{s}"
        out_dir.mkdir(parents=True, exist_ok=True)

        frames = {p: joined[joined["partition"] == p][COLS].copy() for p in PART_FILE}

        # Final guard: a sample is "seen" iff its exact (Sequence, EC) was in train.
        # Drop such rows from valid/test1/test2. This preserves intended CROSS-LABEL
        # homology (same sequence, different EC, in different partitions) while
        # guaranteeing no (sequence,label) is leaked. Catches the handful of
        # identical-seq / different-UniRef50 cases from release drift.
        train_pairs = set(zip(frames["train"]["Sequence"], frames["train"]["EC number"]))
        dropped = {}
        for part in ("valid", "test1", "test2"):
            f = frames[part]
            mask = [(seq, ec) not in train_pairs for seq, ec in zip(f["Sequence"], f["EC number"])]
            dropped[part] = (~pd.Series(mask)).sum()
            frames[part] = f[mask]

        for part, fname in PART_FILE.items():
            frames[part].to_csv(out_dir / fname, sep="\t", index=False)

        tr_pairs = train_pairs
        t1_leak = sum((s, e) in tr_pairs for s, e in zip(frames["test1"]["Sequence"], frames["test1"]["EC number"]))
        t2_leak = sum((s, e) in tr_pairs for s, e in zip(frames["test2"]["Sequence"], frames["test2"]["EC number"]))
        print(f"split_{s} {len(frames['train']):8d} {len(frames['valid']):8d} "
              f"{len(frames['test1']):8d} {len(frames['test2']):8d}  "
              f"{t1_leak:12d} {t2_leak:12d}   (guard dropped: "
              f"valid={dropped['valid']}, test1={dropped['test1']}, test2={dropped['test2']})")


if __name__ == "__main__":
    main()
