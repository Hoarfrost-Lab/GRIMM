# HoarfrostLab Datasets Index

A curated collection of datasets published under **HoarfrostLab**.

| Dataset | Link |
|--------|------|
| `GRIMM` | [link](https://huggingface.co/datasets/HoarfrostLab/GRIMM) |

## GRIMM-EC versions

The HuggingFace repo hosts two versions of the EC splits:

- **`EC_v2/`** — the corrected, leakage-free release. **Use this for all new work.**
  Per-protein SwissProt sequences (UniProt release 2025_02), low-support labels split
  by UniRef50 cluster, independent shuffled folds, true open-set `test2`.
- **`EC_v1/`** — the original release (used by the v1 preprint and parallel works),
  retained for reproducibility. Its `test2` has been corrected to be fully disjoint
  from `train`/`validation`/`test1`; `train`/`validation`/`test1` are unchanged.
  See the dataset card for known limitations. Prefer v2.

Both versions provide `amino_acids/` and `nucleotides/`, each with
`split_1`…`split_5` and `{train, validation, test1, test2}`. Files are tab-separated:
v2 uses `.tsv`; v1 keeps the original `.csv` extension (tab-separated). The custom split
is not a true k-fold but preserves UniRef50 clusters; use the 5 folds as independent
models (individually or as an ensemble).

## Recreate GRIMM-EC v2

    # 1. Sequences + cross-references from UniProt release 2025_02
    #    Download uniprot_sprot-only2025_02.tar.gz from
    #    ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2025_02/knowledgebase/
    #    and extract uniprot_sprot.fasta.gz and uniprot_sprot.dat.gz
    #
    # 2. Build the base table (reuses existing UniRef50 IDs from the v1 release,
    #    joins real per-protein sequences, extracts EMBL CDS ids, normalizes EC labels)
    python code/scripts/combine_data.py \
        --v1-aa-root <EC_v1/amino_acids> --fasta-gz <...fasta.gz> --dat-gz <...dat.gz>
    #
    # 3. Amino-acid splits (by-cluster, seeded, independent folds)
    python code/scripts/create_data_splits_custom.py
    #
    # 4. Nucleotide splits (reformats the AA splits to CDS)
    python code/scripts/create_nuc_splits.py

## Patch GRIMM-EC v1 test-2 (decontamination)

    # Removes from each fold's test2 any row whose accession or sequence appears in
    # train/validation/test1, yielding a fully-disjoint open-set test2.
    python code/scripts/decontaminate_test2.py --in-root <EC_v1> --out-root <patched>

> Note: the original `combine_data.py` / `create_data_splits_custom.py` /
> `create_nuc_splits.py` (which produced v1 and contained the split-leakage bugs)
> are preserved in this repository's git history.
