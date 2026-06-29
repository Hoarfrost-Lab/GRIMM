"""GRIMM-EC nucleotide splits (v2).

Reformats the amino-acid splits into nucleotide (CDS) splits, following the same
accession-level partitioning so the nucleotide splits inherit the clean, by-cluster,
leakage-free structure. Each accession is expanded to its EMBL CDS protein ids
(release-2025_02) and the CDS nucleotide sequence is reused from the v1 nucleotide
release (same accessions -> high coverage).

Fixes vs the original create_nuc_splits.py (in git history): a fresh frame is built
per partition (the original accumulated rows across files, making each output a
cumulative superset), CDS are de-duplicated across the split (train-first), and a
(Sequence, EC) guard drops any eval CDS whose sequence+label is in train.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

PART_FILE = {"train": "train.tsv", "valid": "validation.tsv",
             "test1": "test1.tsv", "test2": "test2.tsv"}
ORDER = ["train", "valid", "test1", "test2"]  # dedup precedence


def load_cds_seqs(v1_nuc_root: Path) -> dict:
    seqs = {}
    for s in range(1, 6):
        for fn in ["train.csv", "valid.csv", "test1.csv", "test2.csv"]:
            f = v1_nuc_root / f"split_{s}" / fn
            if f.exists():
                d = pd.read_csv(f, sep="\t")
                for e, seq in zip(d["Entry"], d["Sequence"]):
                    seqs.setdefault(e, seq)
    return seqs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--aa-root", type=Path, default=Path("v2_build/data/GRIMM_v2/amino_acids"))
    p.add_argument("--cds-map", type=Path, default=Path("v2_build/data/accession_to_cds.tsv"))
    p.add_argument("--v1-nuc-root", type=Path,
                   default=Path("/Users/adrienne/Projects/TEACUP/data/external_benchmarks/"
                                "GRIMM-base/hf_dataset/EC/nucleotides"))
    p.add_argument("--out-root", type=Path, default=Path("v2_build/data/GRIMM_v2/nucleotides"))
    args = p.parse_args()

    cmap = pd.read_csv(args.cds_map, sep="\t").fillna("")
    acc2cds = {e: c.split(";") for e, c in zip(cmap["Entry"], cmap["EMBL_CDS"]) if c}
    cds_seq = load_cds_seqs(args.v1_nuc_root)
    print(f"accession->CDS: {len(acc2cds)} | CDS->seq: {len(cds_seq)}")

    print(f"\n{'split':7} {'train':>8} {'valid':>8} {'test1':>8} {'test2':>8}  {'t1_leak':>8} {'t2_leak':>8}")
    for s in range(1, 6):
        # read AA splits -> accession partition + EC
        part_ec = {}
        for part, fname in PART_FILE.items():
            d = pd.read_csv(args.aa_root / f"split_{s}" / fname, sep="\t")
            for e, ec in zip(d["Entry"], d["EC number"]):
                part_ec[e] = (part, ec)

        seen, rows = set(), {p: [] for p in PART_FILE}
        for part in ORDER:
            for acc, (pp, ec) in part_ec.items():
                if pp != part:
                    continue
                for cid in acc2cds.get(acc, []):
                    if cid in seen or cid not in cds_seq:
                        continue
                    seen.add(cid)
                    rows[part].append((cid, ec, cds_seq[cid]))

        frames = {p: pd.DataFrame(rows[p], columns=["Entry", "EC number", "Sequence"]) for p in PART_FILE}
        train_pairs = set(zip(frames["train"]["Sequence"], frames["train"]["EC number"]))
        for part in ("valid", "test1", "test2"):
            f = frames[part]
            frames[part] = f[[(x, y) not in train_pairs for x, y in zip(f["Sequence"], f["EC number"])]]

        out_dir = args.out_root / f"split_{s}"
        out_dir.mkdir(parents=True, exist_ok=True)
        for part, fname in PART_FILE.items():
            frames[part].to_csv(out_dir / fname, sep="\t", index=False)

        t1 = sum((x, y) in train_pairs for x, y in zip(frames["test1"]["Sequence"], frames["test1"]["EC number"]))
        t2 = sum((x, y) in train_pairs for x, y in zip(frames["test2"]["Sequence"], frames["test2"]["EC number"]))
        print(f"split_{s} {len(frames['train']):8d} {len(frames['valid']):8d} "
              f"{len(frames['test1']):8d} {len(frames['test2']):8d}  {t1:8d} {t2:8d}")


if __name__ == "__main__":
    main()
