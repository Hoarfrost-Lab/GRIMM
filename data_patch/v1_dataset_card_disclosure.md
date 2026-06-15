<!-- Disclosure block to merge into the HoarfrostLab/GRIMM dataset card (v1).
     Draft for PI review — not yet uploaded. -->

## ⚠️ Version note & known issues (please read before use)

This repository hosts **GRIMM v1** (the splits used in the original preprint and in
downstream work). A corrected **v2** is available at `EC_v2/` (see below). We recommend
**v2 for all new work**. v1 is retained, with the correction below, for reproducibility
of already-published results.

### Correction applied to v1 (test-2)

The originally released `test2` ("open-set", out-of-distribution) split inadvertently
included sequences that were also present in `train`: ~80% of `test2` rows shared an
accession (and sequence) with `train`. This was a file-writing error in the split
script (the full orphan pool was written to `test2` instead of the held-out subset).

`test2` has been **corrected to be fully disjoint** from `train`, `validation`, and
`test1` (by both accession and sequence). The corrected `test2` is the genuine
out-of-distribution set (labels absent from training). **`train`, `validation`, and
`test1` are unchanged** from the original release, so results computed on those splits
are unaffected.

### Known limitations remaining in v1 train/validation/test1 (fixed in v2)

These are retained in v1 for citation stability and are documented here for transparency:

1. **Amino-acid sequences are UniRef50 cluster representatives, not per-protein
   SwissProt sequences.** Every accession in a UniRef50 cluster carries that cluster's
   representative sequence (so distinct accessions can share an identical AA sequence).
   The preprint describes per-protein SwissProt sequences; the released v1 AA data does
   not match this. The **nucleotide** modality uses real per-CDS sequences and is not
   affected.
2. **~6% of AA `test1` rows share an exact sequence with `train`.** Labels (EC numbers)
   with only 1–2 UniRef50 clusters were split by individual sequence rather than by
   cluster, so some clusters appear in both `train` and `test1`. Combined with (1) this
   surfaces as ~5.8% of AA `test1` rows (and ~4.1% of `validation`) sharing an exact
   sequence with `train`. The nucleotide modality shows only ~0.4%.

### GRIMM v2 (`EC_v2/`)

v2 regenerates the EC splits to match the documented method and removes the issues above:
- real per-protein SwissProt sequences (UniProt release 2025_02);
- low-support labels split **by UniRef50 cluster** (no by-row splitting);
- independent, shuffled 5-fold partitions;
- `test2` held-out only (true open-set);
- verified **0** `(sequence, EC)` overlap between `train` and any evaluation split.
