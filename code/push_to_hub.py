import os
from huggingface_hub import login, HfApi

# -----------------------------
ORG_NAME = "HoarfrostLab"
DATASETS_ROOT = "data/Stratified/custom/"
HF_TOKEN = os.getenv("HF_TOKEN")
DATASET_NAME = "GRIMM"
INDEX_REPO = f"{ORG_NAME}/hoarfrostlab-datasets-index"

# Expected files per type
AMINO_FILES = ["train.csv", "valid.csv", "test1.csv", "test2.csv"]
NUCLEO_FILES = ["train_reformated.csv", "valid_reformated.csv", "test1_reformated.csv", "test2_reformated.csv"]

# -----------------------------
login(token=HF_TOKEN)
api = HfApi()

# Create the dataset repo if it doesn't exist
api.create_repo(
    repo_id=f"{ORG_NAME}/{DATASET_NAME}",
    repo_type="dataset",
    exist_ok=True,
    token=HF_TOKEN
)

# -----------------------------
# Iterate over folders and upload CSVs
for folder in os.listdir(DATASETS_ROOT):
    folder_path = os.path.join(DATASETS_ROOT, folder)
    if not os.path.isdir(folder_path):
        continue

    # Determine type and expected files
    if folder.startswith("nuc_"):
        dtype = "nucleotides"
        EXPECTED_FILES = NUCLEO_FILES
        split_name = folder.replace("nuc_", "")
    elif folder.startswith("split_"):
        dtype = "amino_acids"
        EXPECTED_FILES = AMINO_FILES
        split_name = folder
    else:
        continue

    # Check for missing files
    files_in_folder = os.listdir(folder_path)
    missing_files = [f for f in EXPECTED_FILES if f not in files_in_folder]
    if missing_files:
        print(f"⚠ Skipping {folder_path}, missing files: {missing_files}")
        continue

    # Upload CSVs
    for f in EXPECTED_FILES:
        subset_name = os.path.splitext(f)[0].lower()
        if subset_name == "valid":
            subset_name = "validation"
        if "_reformated" in subset_name:
            subset_name = subset_name.replace("_reformated", "")

        # Use the original CSV path directly
        csv_path = os.path.join(folder_path, f)

        # Construct HF path
        repo_path = f"EC/{dtype}/{split_name}/{subset_name}.csv"
        api.upload_file(
            path_or_fileobj=csv_path,
            path_in_repo=repo_path,
            repo_id=f"{ORG_NAME}/{DATASET_NAME}",
            repo_type="dataset",
            token=HF_TOKEN
        )
        print(f"✅ Uploaded {repo_path}")

# -----------------------------
# Update index
index_md = (
    "# HoarfrostLab Datasets Index\n\n"
    "A curated collection of datasets published under **HoarfrostLab**.\n\n"
    "| Dataset | Link |\n"
    "|--------|------|\n"
)
index_md += f"| `{DATASET_NAME}` | [link](https://huggingface.co/datasets/{ORG_NAME}/{DATASET_NAME}) |\n"

api.create_repo(
    repo_id=INDEX_REPO,
    repo_type="dataset",
    exist_ok=True,
    token=HF_TOKEN
)
with open("README.md", "w", encoding="utf-8") as f:
    f.write(index_md)

api.upload_file(
    path_or_fileobj="README.md",
    path_in_repo="README.md",
    repo_id=INDEX_REPO,
    repo_type="dataset",
    token=HF_TOKEN,
)

print("✅ Nested dataset pushed and index updated:")
print(f"https://huggingface.co/datasets/{ORG_NAME}/{DATASET_NAME}")

