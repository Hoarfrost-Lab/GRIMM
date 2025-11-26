# HoarfrostLab Datasets Index

A curated collection of datasets published under **HoarfrostLab**.

| Dataset | Link |
|--------|------|
| `GRIMM` | [link](https://huggingface.co/datasets/HoarfrostLab/GRIMM) |

To recreate Stratified data:

    Download links in links.txt
    Run code/scripts/id_mapping.py to get chunks for Swissprot/Uniref50 id mapping
    Run code/scripts/combine_data.py to get the swissprot_full.tsv
    Run appropriate "create_data_splits_{}.py" to get the appropriate data splits (either biological pseudo-kfold (custom) or scikit k-fold (scikit). The custom version is not truly a k-fold but preserves uniref50 clusters. It must be used such that 5 independent models are trained and infered individually or via ensemble.

For nucleotides:

    Run id_mapping.py again for EMBL_CDS ids
    Run add_nuc_to_swissprot.py to combine files into new swissprot
    Run download.py to grab the sequences from ENA and pipe to fasta (courtesy: Soumya)
    Run create_nuc_splits.py to reformat ids in split_[1-5] into nuc_split_[1-5]

