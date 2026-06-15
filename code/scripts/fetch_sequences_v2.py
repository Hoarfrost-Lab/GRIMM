"""Fetch real per-protein SwissProt amino-acid sequences for the v2 base table (Bug C fix).

v1 stored UniRef50 *representative* sequences (combine_data.py merged on UniRef50).
v2 uses each accession's own SwissProt sequence, matching the preprint.

Downloads the SwissProt FASTA once, builds accession->sequence, maps the base table,
and REST-fetches any misses. Records coverage so gaps are explicit (no silent drops).
"""
from __future__ import annotations

import argparse
import gzip
import ssl
import urllib.request
from pathlib import Path

import certifi
import pandas as pd

CTX = ssl.create_default_context(cafile=certifi.where())
FASTA_URL = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
             "knowledgebase/complete/uniprot_sprot.fasta.gz")
REST = "https://rest.uniprot.org/uniprotkb/{}.fasta"


def download(url: str, dest: Path) -> None:
    if dest.exists():
        print(f"already have {dest} ({dest.stat().st_size/1e6:.1f} MB)")
        return
    print(f"downloading {url} -> {dest}")
    with urllib.request.urlopen(url, context=CTX, timeout=600) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    print(f"done ({dest.stat().st_size/1e6:.1f} MB)")


def parse_fasta_gz(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    acc = None
    buf: list[str] = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf)
                # header: >sp|ACCESSION|NAME ...
                acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if acc is not None:
            seqs[acc] = "".join(buf)
    return seqs


def rest_fetch(acc: str) -> str | None:
    try:
        with urllib.request.urlopen(REST.format(acc), context=CTX, timeout=30) as r:
            txt = r.read().decode()
        return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">")) or None
    except Exception:
        return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base-table", type=Path, default=Path("v2_build/data/base_table.tsv"))
    p.add_argument("--fasta", type=Path, default=Path("v2_build/data/uniprot_sprot.fasta.gz"))
    p.add_argument("--out", type=Path, default=Path("v2_build/data/base_table_with_seq.tsv"))
    p.add_argument("--max-rest", type=int, default=5000, help="cap REST fallback fetches")
    args = p.parse_args()

    download(FASTA_URL, args.fasta)
    print("parsing FASTA ...")
    seqs = parse_fasta_gz(args.fasta)
    print(f"  {len(seqs)} SwissProt sequences parsed")

    df = pd.read_csv(args.base_table, sep="\t")
    df["Sequence"] = df["Entry"].map(seqs)
    missing = df[df["Sequence"].isna()]["Entry"].tolist()
    print(f"mapped {len(df)-len(missing)}/{len(df)} from FASTA; {len(missing)} misses")

    if missing:
        n = min(len(missing), args.max_rest)
        print(f"REST fallback for {n} misses (cap {args.max_rest}) ...")
        fetched = {a: rest_fetch(a) for a in missing[:n]}
        df.loc[df["Entry"].isin(fetched), "Sequence"] = df["Entry"].map(fetched)
        still = df["Sequence"].isna().sum()
        if len(missing) > args.max_rest:
            print(f"  NOTE: {len(missing)-args.max_rest} misses left unfetched (cap)")
        print(f"  still missing after REST: {still}")

    got = df["Sequence"].notna().sum()
    print(f"FINAL coverage: {got}/{len(df)} ({100*got/len(df):.2f}%)")
    df.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
