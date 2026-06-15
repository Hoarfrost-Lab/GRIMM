"""Build v2 nucleotide splits from the verified accession-level v2 partitions.

Maps each accession to its EMBL CDS protein ids (release-2025_02, from
accession_to_cds.tsv), looks up the CDS nucleotide sequence (reused from the v1
nucleotide release — same accessions, so coverage is high), and writes per-partition
CDS rows following the v2 (clean, by-cluster) partition assignment.

Fixes vs v1 create_nuc_splits.py:
  Bug E: fresh frame per partition (no cumulative supersets); dedup CDS across the
         split with train-first precedence.
  Bug B/D/A: inherited from the clean v2 accession splits.
Plus a (Sequence, EC) guard so no CDS sequence+label is shared train<->eval.

Inputs:
  v2_build/data/v2_splits/split_*/assignments.tsv   (Entry, EC number, UniRef50, partition)
  v2_build/data/accession_to_cds.tsv                (Entry -> EMBL_CDS ids)
  v1 nucleotide cache (union)                        (CDS id -> nucleotide sequence)
Output:
  v2_build/data/GRIMM_v2/nucleotides/split_*/{train,validation,test1,test2}.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

PART_FILE = {"train": "train.csv", "valid": "validation.csv",
             "test1": "test1.csv", "test2": "test2.csv"}
ORDER = ["train", "valid", "test1", "test2"]  # dedup precedence (train first)


def load_cds_seqs(v1_nuc_root: Path) -> dict:
    """CDS id -> nucleotide sequence, unioned over all v1 nucleotide partitions."""
    seqs: dict = {}
    for s in range(1, 6):
        for fn in ["train.csv", "valid.csv", "test1.csv", "test2.csv"]:
            f = v1_nuc_root / f"split_{s}" / fn
            if not f.exists():
                continue
            d = pd.read_csv(f, sep="\t")
            for entry, seq in zip(d["Entry"], d["Sequence"]):
                seqs.setdefault(entry, seq)
    return seqs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--splits-root", type=Path, default=Path("v2_build/data/v2_splits"))
    p.add_argument("--cds-map", type=Path, default=Path("v2_build/data/accession_to_cds.tsv"))
    p.add_argument("--v1-nuc-root", type=Path,
                   default=Path("/Users/adrienne/Projects/TEACUP/data/external_benchmarks/"
                                "GRIMM-base/hf_dataset/EC/nucleotides"))
    p.add_argument("--out-root", type=Path, default=Path("v2_build/data/GRIMM_v2/nucleotides"))
    p.add_argument("--splits", type=int, default=5)
    args = p.parse_args()

    cds_map = pd.read_csv(args.cds_map, sep="\t").fillna("")
    acc2cds = {e: c.split(";") for e, c in zip(cds_map["Entry"], cds_map["EMBL_CDS"]) if c}
    print(f"accession->CDS map: {len(acc2cds)} accessions")
    cds_seq = load_cds_seqs(args.v1_nuc_root)
    print(f"CDS id->sequence (from v1 nuc): {len(cds_seq)} CDS")

    print(f"\n{'split':7} {'train':>8} {'valid':>8} {'test1':>8} {'test2':>8}  {'t1_leak':>8} {'t2_leak':>8}")
    for s in range(1, args.splits + 1):
        a = pd.read_csv(args.splits_root / f"split_{s}" / "assignments.tsv", sep="\t")
        part_of = dict(zip(a["Entry"], a["partition"]))
        ec_of = dict(zip(a["Entry"], a["EC number"]))

        seen_cds: set = set()
        rows = {p: [] for p in PART_FILE}
        for part in ORDER:  # train first so a shared CDS resolves to train
            accs = a[a["partition"] == part]["Entry"]
            for acc in accs:
                for cid in acc2cds.get(acc, []):
                    if cid in seen_cds or cid not in cds_seq:
                        continue
                    seen_cds.add(cid)
                    rows[part].append((cid, ec_of[acc], cds_seq[cid]))

        frames = {p: pd.DataFrame(rows[p], columns=["Entry", "EC number", "Sequence"])
                  for p in PART_FILE}
        # (Sequence, EC) guard
        train_pairs = set(zip(frames["train"]["Sequence"], frames["train"]["EC number"]))
        for part in ("valid", "test1", "test2"):
            f = frames[part]
            mask = [(seq, ec) not in train_pairs for seq, ec in zip(f["Sequence"], f["EC number"])]
            frames[part] = f[mask]

        out_dir = args.out_root / f"split_{s}"
        out_dir.mkdir(parents=True, exist_ok=True)
        for part, fname in PART_FILE.items():
            frames[part].to_csv(out_dir / fname, sep="\t", index=False)

        t1_leak = sum((x, y) in train_pairs for x, y in zip(frames["test1"]["Sequence"], frames["test1"]["EC number"]))
        t2_leak = sum((x, y) in train_pairs for x, y in zip(frames["test2"]["Sequence"], frames["test2"]["EC number"]))
        print(f"split_{s} {len(frames['train']):8d} {len(frames['valid']):8d} "
              f"{len(frames['test1']):8d} {len(frames['test2']):8d}  {t1_leak:8d} {t2_leak:8d}")


if __name__ == "__main__":
    main()
