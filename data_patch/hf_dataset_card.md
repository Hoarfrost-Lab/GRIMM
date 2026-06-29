---
license: cc-by-4.0
language:
  - en
tags:
  - biology
  - proteins
  - enzymes
  - ec-number
  - genomics
  - protein-function-prediction
pretty_name: GRIMM (EC) — Genomic Representation Inference for Microbial Metabolism
task_categories:
  - text-classification
size_categories:
  - 100K<n<1M
configs:
  - config_name: EC_v2_amino_acids
    default: true
    data_files:
      - split: split1_train
        path: EC_v2/amino_acids/split_1/train.tsv
      - split: split1_validation
        path: EC_v2/amino_acids/split_1/validation.tsv
      - split: split1_test1
        path: EC_v2/amino_acids/split_1/test1.tsv
      - split: split1_test2
        path: EC_v2/amino_acids/split_1/test2.tsv
      - split: split2_train
        path: EC_v2/amino_acids/split_2/train.tsv
      - split: split2_validation
        path: EC_v2/amino_acids/split_2/validation.tsv
      - split: split2_test1
        path: EC_v2/amino_acids/split_2/test1.tsv
      - split: split2_test2
        path: EC_v2/amino_acids/split_2/test2.tsv
      - split: split3_train
        path: EC_v2/amino_acids/split_3/train.tsv
      - split: split3_validation
        path: EC_v2/amino_acids/split_3/validation.tsv
      - split: split3_test1
        path: EC_v2/amino_acids/split_3/test1.tsv
      - split: split3_test2
        path: EC_v2/amino_acids/split_3/test2.tsv
      - split: split4_train
        path: EC_v2/amino_acids/split_4/train.tsv
      - split: split4_validation
        path: EC_v2/amino_acids/split_4/validation.tsv
      - split: split4_test1
        path: EC_v2/amino_acids/split_4/test1.tsv
      - split: split4_test2
        path: EC_v2/amino_acids/split_4/test2.tsv
      - split: split5_train
        path: EC_v2/amino_acids/split_5/train.tsv
      - split: split5_validation
        path: EC_v2/amino_acids/split_5/validation.tsv
      - split: split5_test1
        path: EC_v2/amino_acids/split_5/test1.tsv
      - split: split5_test2
        path: EC_v2/amino_acids/split_5/test2.tsv
  - config_name: EC_v2_nucleotides
    data_files:
      - split: split1_train
        path: EC_v2/nucleotides/split_1/train.tsv
      - split: split1_validation
        path: EC_v2/nucleotides/split_1/validation.tsv
      - split: split1_test1
        path: EC_v2/nucleotides/split_1/test1.tsv
      - split: split1_test2
        path: EC_v2/nucleotides/split_1/test2.tsv
      - split: split2_train
        path: EC_v2/nucleotides/split_2/train.tsv
      - split: split2_validation
        path: EC_v2/nucleotides/split_2/validation.tsv
      - split: split2_test1
        path: EC_v2/nucleotides/split_2/test1.tsv
      - split: split2_test2
        path: EC_v2/nucleotides/split_2/test2.tsv
      - split: split3_train
        path: EC_v2/nucleotides/split_3/train.tsv
      - split: split3_validation
        path: EC_v2/nucleotides/split_3/validation.tsv
      - split: split3_test1
        path: EC_v2/nucleotides/split_3/test1.tsv
      - split: split3_test2
        path: EC_v2/nucleotides/split_3/test2.tsv
      - split: split4_train
        path: EC_v2/nucleotides/split_4/train.tsv
      - split: split4_validation
        path: EC_v2/nucleotides/split_4/validation.tsv
      - split: split4_test1
        path: EC_v2/nucleotides/split_4/test1.tsv
      - split: split4_test2
        path: EC_v2/nucleotides/split_4/test2.tsv
      - split: split5_train
        path: EC_v2/nucleotides/split_5/train.tsv
      - split: split5_validation
        path: EC_v2/nucleotides/split_5/validation.tsv
      - split: split5_test1
        path: EC_v2/nucleotides/split_5/test1.tsv
      - split: split5_test2
        path: EC_v2/nucleotides/split_5/test2.tsv
  - config_name: EC_v1_amino_acids
    data_files:
      - split: split1_train
        path: EC_v1/amino_acids/split_1/train.csv
      - split: split1_validation
        path: EC_v1/amino_acids/split_1/validation.csv
      - split: split1_test1
        path: EC_v1/amino_acids/split_1/test1.csv
      - split: split1_test2
        path: EC_v1/amino_acids/split_1/test2.csv
      - split: split2_train
        path: EC_v1/amino_acids/split_2/train.csv
      - split: split2_validation
        path: EC_v1/amino_acids/split_2/validation.csv
      - split: split2_test1
        path: EC_v1/amino_acids/split_2/test1.csv
      - split: split2_test2
        path: EC_v1/amino_acids/split_2/test2.csv
      - split: split3_train
        path: EC_v1/amino_acids/split_3/train.csv
      - split: split3_validation
        path: EC_v1/amino_acids/split_3/validation.csv
      - split: split3_test1
        path: EC_v1/amino_acids/split_3/test1.csv
      - split: split3_test2
        path: EC_v1/amino_acids/split_3/test2.csv
      - split: split4_train
        path: EC_v1/amino_acids/split_4/train.csv
      - split: split4_validation
        path: EC_v1/amino_acids/split_4/validation.csv
      - split: split4_test1
        path: EC_v1/amino_acids/split_4/test1.csv
      - split: split4_test2
        path: EC_v1/amino_acids/split_4/test2.csv
      - split: split5_train
        path: EC_v1/amino_acids/split_5/train.csv
      - split: split5_validation
        path: EC_v1/amino_acids/split_5/validation.csv
      - split: split5_test1
        path: EC_v1/amino_acids/split_5/test1.csv
      - split: split5_test2
        path: EC_v1/amino_acids/split_5/test2.csv
  - config_name: EC_v1_nucleotides
    data_files:
      - split: split1_train
        path: EC_v1/nucleotides/split_1/train.csv
      - split: split1_validation
        path: EC_v1/nucleotides/split_1/validation.csv
      - split: split1_test1
        path: EC_v1/nucleotides/split_1/test1.csv
      - split: split1_test2
        path: EC_v1/nucleotides/split_1/test2.csv
      - split: split2_train
        path: EC_v1/nucleotides/split_2/train.csv
      - split: split2_validation
        path: EC_v1/nucleotides/split_2/validation.csv
      - split: split2_test1
        path: EC_v1/nucleotides/split_2/test1.csv
      - split: split2_test2
        path: EC_v1/nucleotides/split_2/test2.csv
      - split: split3_train
        path: EC_v1/nucleotides/split_3/train.csv
      - split: split3_validation
        path: EC_v1/nucleotides/split_3/validation.csv
      - split: split3_test1
        path: EC_v1/nucleotides/split_3/test1.csv
      - split: split3_test2
        path: EC_v1/nucleotides/split_3/test2.csv
      - split: split4_train
        path: EC_v1/nucleotides/split_4/train.csv
      - split: split4_validation
        path: EC_v1/nucleotides/split_4/validation.csv
      - split: split4_test1
        path: EC_v1/nucleotides/split_4/test1.csv
      - split: split4_test2
        path: EC_v1/nucleotides/split_4/test2.csv
      - split: split5_train
        path: EC_v1/nucleotides/split_5/train.csv
      - split: split5_validation
        path: EC_v1/nucleotides/split_5/validation.csv
      - split: split5_test1
        path: EC_v1/nucleotides/split_5/test1.csv
      - split: split5_test2
        path: EC_v1/nucleotides/split_5/test2.csv
---

# GRIMM-EC

GRIMM is a benchmark for predicting enzyme function (EC number) from biological
sequence. Sequences are reviewed (SwissProt) prokaryotic proteins with EC
annotations; partitions are stratified per label by UniRef50 cluster so that
homologous sequences do not leak between train and evaluation splits.

Two modalities are provided: **amino acids** (per-protein SwissProt sequences) and
**nucleotides** (per-CDS sequences from ENA).

## ⭐ Which version to use

| Path | Status | Use it? |
|------|--------|---------|
| **`EC_v2/`** | corrected, leakage-free, matches the current preprint | ✅ **Yes — all new work** |
| `EC_v1/`     | original release (v1 preprint + parallel works); legacy | only to reproduce already-published v1 results |

## Structure

```
EC_v2/   amino_acids/   split_1 … split_5 / {train, validation, test1, test2}.tsv
         nucleotides/   split_1 … split_5 / {train, validation, test1, test2}.tsv
EC_v1/   amino_acids/   split_1 … split_5 / {train, validation, test1, test2}.csv
         nucleotides/   split_1 … split_5 / {train, validation, test1, test2}.csv
```

All files are **tab-separated**. v2 uses the `.tsv` extension; v1 retains the original
`.csv` extension (tab-separated despite the name) for release stability — load v1 with
`pd.read_csv(path, sep="\t")`.

Columns — amino acids: `Entry, EC number, EMBL, RefSeq, UniRef50, UniRef90, UniRef100, Sequence`;
nucleotides: `Entry` (EMBL CDS id), `EC number`, `Sequence`.

The 5 folds are **not** a standard k-fold: each is an independent train/valid/test
partition that preserves UniRef50 clusters. Train and evaluate the 5 folds as
independent models (individually or as an ensemble), not as rotating CV folds.

### Splits

- **train / validation / test1** — closed-set: evaluation sequences whose labels also
  appear in training, but held out by UniRef50 cluster.
- **test2** — open-set: sequences from labels **absent from training** (out-of-distribution).

### Labels

EC numbers (4th level). Proteins with multiple EC annotations are kept as a single
**compound label string** (e.g. `1.1.99.1;1.2.1.8`), distinct from its component
labels — they are **not** expanded into separate rows.

### Sizes (GRIMM-EC v2, average per fold)

| Modality | train | validation | test1 | test2 |
|----------|------:|-----------:|------:|------:|
| amino acids | ~178,053 | ~28,719 | ~29,689 | ~959 |
| nucleotides | ~251,745 | ~42,557 | ~45,185 | ~1,755 |

237,421 proteins · 6,393 EC labels (1,321 compound) · 65,996 UniRef50 clusters.
Sequences from UniProt release **2025_02**.

## How GRIMM-EC v2 is built (and how it differs from v1)

v2 reuses v1's UniRef50 cluster assignments but regenerates the splits to match the
documented method:

- **Per-protein SwissProt sequences** (release 2025_02) — v1's AA data instead held
  the UniRef50 *representative* sequence.
- **Low-support labels split by UniRef50 cluster** — labels with 1–2 clusters are
  partitioned by whole cluster (2 clusters → 1 train / 1 test1; 1 cluster → orphan,
  ~80% train / ~20% test2 across folds), not by individual sequence as in v1.
- **Independent, shuffled folds**; **seeded** for reproducibility.
- **`test2` is held-out only** (true open-set) — v1 inadvertently also wrote the
  held-out orphans into train.
- EC labels normalized (stray whitespace stripped).

**Verified for v2 (all 5 folds, both modalities):** 0 `(sequence, EC)` overlap between
`train` and any evaluation split; 0 accession overlap between splits; `test2` labels
absent from train. Identical sequences carrying *different* EC labels may appear in
different splits — this is intended cross-label difficulty under per-label
stratification, not leakage.

## Known limitations of GRIMM-EC v1 (fixed in v2)

`EC_v1/` is retained for reproducibility of already-published results. Its `test2`
has been **corrected** to be fully disjoint from `train`/`validation`/`test1`;
`train`/`validation`/`test1` are **unchanged** from the original release. Remaining
v1 limitations (all fixed in v2):

1. AA sequences are UniRef50 **representatives**, not per-protein SwissProt sequences
   (the nucleotide modality used real per-CDS sequences and is unaffected).
2. ~5.8% of AA `test1` rows (and ~4.1% of `validation`) share an exact sequence with
   `train`, because 1–2 cluster labels were split by sequence rather than by cluster
   (nucleotides: ~0.4%).

See the repository / preprint for full details.

## Citation

> Hoarfrost et al. GRIMM: Genomic Representation Inference for Microbial Metabolism. (preprint)

Code: https://github.com/Hoarfrost-Lab/grimm
