"""Extract May-2025 (release-2025_02) SwissProt sequences + EMBL-CDS ids from the
release tarball, restricted to the v2 base-table accessions.

Streams the two needed members out of the 1.7 GB tar.gz without extracting to disk:
  - uniprot_sprot.fasta.gz  -> accession -> real AA sequence (Bug C, release-pinned)
  - uniprot_sprot.dat.gz     -> accession -> EMBL CDS protein ids (for nucleotides)

Outputs:
  v2_build/data/base_table_with_seq.tsv   (base table + release-2025_02 Sequence)
  v2_build/data/accession_to_cds.tsv      (Entry -> semicolon EMBL CDS protein ids)
"""
from __future__ import annotations

import argparse
import gzip
import io
import tarfile
from pathlib import Path

import pandas as pd


def find_member(tar: tarfile.TarFile, suffix: str) -> tarfile.TarInfo:
    for m in tar.getmembers():
        if m.name.endswith(suffix):
            return m
    raise KeyError(suffix)


def parse_fasta(stream, wanted: set) -> dict:
    seqs, acc, buf = {}, None, []
    gz = stream
    for raw in io.TextIOWrapper(gz, encoding="ascii"):
        if raw.startswith(">"):
            if acc is not None and acc in wanted:
                seqs[acc] = "".join(buf)
            acc = raw[1:].split("|")[1] if "|" in raw else raw[1:].split()[0]
            buf = []
        else:
            buf.append(raw.strip())
    if acc is not None and acc in wanted:
        seqs[acc] = "".join(buf)
    return seqs


def parse_dat_cds(stream, wanted: set) -> dict:
    """accession -> list of EMBL CDS protein ids (3rd field of DR EMBL lines)."""
    cds: dict = {}
    gz = stream
    acc = None
    ids: list = []
    for raw in io.TextIOWrapper(gz, encoding="ascii", errors="replace"):
        if raw.startswith("AC "):
            if acc is None:  # first AC line of the entry = primary accession
                acc = raw[5:].split(";")[0].strip()
        elif raw.startswith("DR   EMBL;"):
            parts = [p.strip() for p in raw.split(";")]
            # DR   EMBL; <nucleotide>; <protein/CDS id>; <moltype>; <status>.
            if len(parts) >= 3 and parts[2] not in ("-", ""):
                ids.append(parts[2])
        elif raw.startswith("//"):
            if acc is not None and acc in wanted and ids:
                cds[acc] = ";".join(dict.fromkeys(ids))  # dedup, keep order
            acc, ids = None, []
    return cds


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta-gz", type=Path, default=Path("v2_build/data/uniprot_sprot.fasta.gz"))
    p.add_argument("--dat-gz", type=Path, default=Path("v2_build/data/uniprot_sprot.dat.gz"))
    p.add_argument("--base-table", type=Path, default=Path("v2_build/data/base_table.tsv"))
    p.add_argument("--out-seq", type=Path, default=Path("v2_build/data/base_table_with_seq.tsv"))
    p.add_argument("--out-cds", type=Path, default=Path("v2_build/data/accession_to_cds.tsv"))
    args = p.parse_args()

    base = pd.read_csv(args.base_table, sep="\t")
    wanted = set(base["Entry"])
    print(f"base table: {len(wanted)} accessions")

    print(f"parsing {args.fasta_gz} ...")
    with gzip.open(args.fasta_gz, "rb") as fh:
        seqs = parse_fasta(fh, wanted)
    print(f"  got {len(seqs)} sequences (of {len(wanted)})")

    print(f"parsing {args.dat_gz} for EMBL CDS ids ...")
    with gzip.open(args.dat_gz, "rb") as fh:
        cds = parse_dat_cds(fh, wanted)
    print(f"  got CDS ids for {len(cds)} accessions")

    base["Sequence"] = base["Entry"].map(seqs)
    cov = base["Sequence"].notna().sum()
    print(f"AA sequence coverage: {cov}/{len(base)} ({100*cov/len(base):.2f}%)")
    base.to_csv(args.out_seq, sep="\t", index=False)
    pd.DataFrame({"Entry": list(cds), "EMBL_CDS": list(cds.values())}).to_csv(
        args.out_cds, sep="\t", index=False)
    print(f"wrote {args.out_seq} and {args.out_cds}")


if __name__ == "__main__":
    main()
