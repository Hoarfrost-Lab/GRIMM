"""Build the GRIMM-EC base table (v2).

v2 reuses the existing UniRef50 cluster assignments from the released v1 dataset
(so cluster IDs are unchanged) but takes each protein's REAL per-protein SwissProt
amino-acid sequence from UniProt release 2025_02 — fixing the v1 issue where the
sequence column held the UniRef50 cluster *representative* (see README / preprint).

It also extracts EMBL CDS protein ids (for the nucleotide splits) from the same
release, and normalizes EC labels (strips stray whitespace in components so the
same EC is not treated as two labels).

Inputs:
  --v1-aa-root   released v1 EC/amino_acids dir (for Entry, EC, EMBL, RefSeq,
                 UniRef50/90/100 metadata + existing UniRef50 IDs)
  --fasta-gz     uniprot_sprot.fasta.gz from release-2025_02 (per-protein AA seqs)
  --dat-gz       uniprot_sprot.dat.gz   from release-2025_02 (EMBL CDS xrefs)

Outputs:
  base_table_with_seq.tsv   Entry, EC number, EMBL, RefSeq, UniRef50/90/100, Sequence
  accession_to_cds.tsv      Entry, EMBL_CDS  (semicolon-separated CDS protein ids)
"""
from __future__ import annotations

import argparse
import gzip
import io
from pathlib import Path

import pandas as pd

META = ["Entry", "EC number", "EMBL", "RefSeq", "UniRef50", "UniRef90", "UniRef100"]


def normalize_ec(label: str) -> str:
    """Strip whitespace from each ';'-separated EC component (e.g. ' 2.3.1.266')."""
    return ";".join(p.strip() for p in str(label).split(";"))


def reconstruct_metadata(v1_aa_root: Path) -> pd.DataFrame:
    frames = []
    for s in range(1, 6):
        for fn in ["train.csv", "validation.csv", "test1.csv", "test2.csv"]:
            f = v1_aa_root / f"split_{s}" / fn
            if f.exists():
                frames.append(pd.read_csv(f, sep="\t")[META])
    full = pd.concat(frames, ignore_index=True).drop_duplicates(subset="Entry")
    full["EC number"] = full["EC number"].map(normalize_ec)
    return full.reset_index(drop=True)


def parse_fasta(path: Path, wanted: set) -> dict:
    seqs, acc, buf = {}, None, []
    with gzip.open(path, "rb") as raw:
        for line in io.TextIOWrapper(raw, encoding="ascii"):
            if line.startswith(">"):
                if acc is not None and acc in wanted:
                    seqs[acc] = "".join(buf)
                acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None and acc in wanted:
        seqs[acc] = "".join(buf)
    return seqs


def parse_dat_cds(path: Path, wanted: set) -> dict:
    cds, acc, ids = {}, None, []
    with gzip.open(path, "rb") as raw:
        for line in io.TextIOWrapper(raw, encoding="ascii", errors="replace"):
            if line.startswith("AC ") and acc is None:
                acc = line[5:].split(";")[0].strip()
            elif line.startswith("DR   EMBL;"):
                parts = [p.strip() for p in line.split(";")]
                if len(parts) >= 3 and parts[2] not in ("-", ""):
                    ids.append(parts[2])
            elif line.startswith("//"):
                if acc in wanted and ids:
                    cds[acc] = ";".join(dict.fromkeys(ids))
                acc, ids = None, []
    return cds


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--v1-aa-root", type=Path,
                   default=Path("/Users/adrienne/Projects/TEACUP/data/external_benchmarks/"
                                "GRIMM-base/hf_dataset/EC/amino_acids"))
    p.add_argument("--fasta-gz", type=Path, default=Path("v2_build/data/uniprot_sprot.fasta.gz"))
    p.add_argument("--dat-gz", type=Path, default=Path("v2_build/data/uniprot_sprot.dat.gz"))
    p.add_argument("--out-seq", type=Path, default=Path("v2_build/data/base_table_with_seq.tsv"))
    p.add_argument("--out-cds", type=Path, default=Path("v2_build/data/accession_to_cds.tsv"))
    args = p.parse_args()

    base = reconstruct_metadata(args.v1_aa_root)
    wanted = set(base["Entry"])
    print(f"base table: {len(base)} accessions, {base['EC number'].nunique()} EC labels "
          f"(normalized), {base['UniRef50'].nunique()} UniRef50 clusters")

    print("parsing release-2025_02 FASTA (real per-protein sequences) ...")
    base["Sequence"] = base["Entry"].map(parse_fasta(args.fasta_gz, wanted))
    cov = base["Sequence"].notna().sum()
    print(f"  AA sequence coverage: {cov}/{len(base)} ({100*cov/len(base):.2f}%)")

    print("parsing release-2025_02 .dat (EMBL CDS ids) ...")
    cds = parse_dat_cds(args.dat_gz, wanted)
    print(f"  CDS ids for {len(cds)} accessions")

    args.out_seq.parent.mkdir(parents=True, exist_ok=True)
    base.to_csv(args.out_seq, sep="\t", index=False)
    pd.DataFrame({"Entry": list(cds), "EMBL_CDS": list(cds.values())}).to_csv(
        args.out_cds, sep="\t", index=False)
    print(f"wrote {args.out_seq} and {args.out_cds}")


if __name__ == "__main__":
    main()
