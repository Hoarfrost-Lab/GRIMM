# GRIMM data-split fix — review & sign-off

**For:** coauthor review before we re-upload to HuggingFace.
**Status:** all changes are on branch `fix/data-split-leakage` (this PR). **Nothing has been pushed to HuggingFace yet** — that step waits for your sign-off.

## TL;DR

We found and fixed a data-split leakage bug in GRIMM-EC. The headline issue: ~80% of the `test2` (open-set) sequences were also present in `train`. Investigating it surfaced a few related issues. We are doing two things:

1. **Patch v1 in place** (minimal): re-upload a corrected `test2` only; **`train`/`validation`/`test1` stay byte-identical**, so our already-submitted papers (which use the AA train/test1 splits) are unaffected. v1 is renamed `EC_v1/` and kept for citation stability.
2. **Publish a clean `EC_v2/`** that fixes everything and matches what the preprint describes. **Use v2 for all new work.**

## The bugs (verified against the live HF dataset)

| | Bug | Impact (measured) |
|--|-----|-------------------|
| A | `test2` written from the full orphan pool instead of the held-out subset | **80%** of `test2` rows also in `train` |
| B | Labels with 1–2 UniRef50 clusters split by *row*, not by *cluster* | **5.8%** of AA `test1` rows share an exact sequence with `train` |
| C | AA `Sequence` column held the UniRef50 *representative*, not the protein's own SwissProt sequence | (amplifies B; preprint says per-protein) |
| D | 5 folds not independent (deterministic cluster ordering) | test1 ~95% identical across folds |
| E/F | nucleotide reformatting bug; random seed declared but never used | reproducibility |

B and C are deviations of the *released data* from the method described in the preprint — i.e. the paper's description is correct; the released data didn't match it. v2 makes the data match the paper.

## What v1 (`EC_v1/`) gets

- `test2` **decontaminated** to be fully disjoint from `train`/`validation`/`test1` (by accession and sequence). 0 residual leakage.
- `train`/`validation`/`test1` **unchanged**.
- Known limitations (B, C) documented in the dataset card; the exact original is still retrievable at HF commit `42b52d6`.

## What v2 (`EC_v2/`) is

Regenerated to match the preprint, **reusing the existing UniRef50 cluster IDs**:
- Real per-protein **SwissProt sequences** (UniProt release **2025_02**), fixing C.
- Low-support labels split **by UniRef50 cluster** (2 clusters → 1 train / 1 test1; 1 cluster → orphan, ~80/20 train/test2 across folds), fixing B.
- **Independent, shuffled, seeded** folds (D, F); `test2` held-out only (A).
- EC labels normalized (stray whitespace stripped — fixed 283 labels that were spuriously split).
- **Multi-EC handling unchanged from v1:** compound annotations (e.g. `1.1.99.1;1.2.1.8`) are kept as a single distinct label, *not* expanded.

### v2 QC (all 5 folds, both modalities)
- **0** `(sequence, EC)` overlap between `train` and any of validation/test1/test2.
- **0** accession overlap between splits.
- **0** `test2` labels present in `train` (true open-set).
- Unique sequences ≫ unique clusters (confirms real per-protein sequences; v1 had them equal).
- Identical sequences with *different* EC labels may appear across splits — intended cross-label difficulty under per-label stratification, not leakage.

### v2 sizes (avg per fold)
| Modality | train | validation | test1 | test2 |
|---|--:|--:|--:|--:|
| amino acids | ~178,053 | ~28,719 | ~29,689 | ~959 |
| nucleotides | ~251,745 | ~42,557 | ~45,185 | ~1,755 |

237,421 proteins · 6,393 EC labels · 65,996 UniRef50 clusters.

## What to look at (links are on this branch)

- **Plan & rationale:** [`grimm_leakage_fix_plan.md`](grimm_leakage_fix_plan.md)
- **v1 dataset-card disclosure (your edits):** [`data_patch/v1_dataset_card_disclosure.md`](data_patch/v1_dataset_card_disclosure.md)
- **Full v2 HuggingFace dataset card** (what will be published): [`data_patch/hf_dataset_card.md`](data_patch/hf_dataset_card.md)
- **Patched v1 `test2` files** (the only v1 change): [`data_patch/v1_test2/`](data_patch/v1_test2/)
- **Pipeline scripts** (under original names; buggy originals are in git history):
  [`code/scripts/combine_data.py`](code/scripts/combine_data.py),
  [`create_data_splits_custom.py`](code/scripts/create_data_splits_custom.py),
  [`create_nuc_splits.py`](code/scripts/create_nuc_splits.py),
  [`decontaminate_test2.py`](code/scripts/decontaminate_test2.py)
- **Deployment script** (dry-run by default): [`code/push_to_hub.py`](code/push_to_hub.py)

> Note: the large v2 CSVs themselves are not in git (they'll appear on HuggingFace under `EC_v2/` after sign-off). The QC numbers and sizes above are from the generated files; happy to share samples directly if useful.

## What needs your sign-off

The HuggingFace upload, which will:
1. Rename the current `EC/` → `EC_v1/` (and overwrite its `test2` with the patched version).
2. Add the clean `EC_v2/` (default config on the Hub).
3. Publish the dataset card.

⚠️ Renaming `EC/` → `EC_v1/` means any code that loads the **`EC/` path from the Hub** will need to point at `EC_v1/` or `EC_v2/`. (Our submitted papers use local copies, so they're unaffected.)

Once you approve, we run `python code/push_to_hub.py --execute`, then merge this branch to `main`.
