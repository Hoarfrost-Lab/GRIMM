"""Deploy the GRIMM dataset correction to HuggingFace (HoarfrostLab/GRIMM).

Workflow:
  1. Rename the live `EC/` tree to `EC_v1/` (server-side copy + delete — no re-upload
     of the large LFS files).
  2. Overwrite `EC_v1/.../test2.csv` with the decontaminated (fully-disjoint) test-2.
  3. Upload the corrected `EC_v2/` tree (amino_acids + nucleotides, 5 folds).
  4. Upload the dataset card (README.md).

SAFETY: dry-run by default — prints the planned operations and changes NOTHING.
Pass --execute to actually commit/push.

  python code/push_to_hub.py                 # dry run (default)
  python code/push_to_hub.py --execute       # perform the deployment

Auth uses the cached HuggingFace token (~/.cache/huggingface/token) or HF_TOKEN.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

from huggingface_hub import (HfApi, hf_hub_download, CommitOperationAdd,
                             CommitOperationCopy, CommitOperationDelete)

REPO = "HoarfrostLab/GRIMM"
RTYPE = "dataset"
MODALITIES = ["amino_acids", "nucleotides"]
FILES = ["train.csv", "validation.csv", "test1.csv", "test2.csv"]


def plan_rename(api, execute):
    # Determine LFS status per file: huggingface_hub 0.16.4 can server-side-COPY
    # only LFS files. Non-LFS files (here, the small test2.csv) are downloaded and
    # re-added instead. Both cases delete the old EC/ path.
    info = api.repo_info(REPO, repo_type=RTYPE, files_metadata=True)
    is_lfs = {s.rfilename: bool(getattr(s, "lfs", None)) for s in info.siblings}
    files = [f for f in is_lfs if f.startswith("EC/")]
    n_lfs = sum(is_lfs[f] for f in files)
    print(f"[1] rename EC/ -> EC_v1/  ({len(files)} files: {n_lfs} LFS copy, "
          f"{len(files)-n_lfs} non-LFS re-add)")
    if not files:
        print("    (no EC/ files found — already renamed?)")
        return
    ops = []
    for f in files:
        dest = "EC_v1/" + f[len("EC/"):]
        if is_lfs[f]:
            ops.append(CommitOperationCopy(src_path_in_repo=f, path_in_repo=dest))
        elif execute:
            local = hf_hub_download(REPO, f, repo_type=RTYPE)
            ops.append(CommitOperationAdd(path_in_repo=dest, path_or_fileobj=local))
        ops.append(CommitOperationDelete(path_in_repo=f))
    if execute:
        retry(lambda: api.create_commit(REPO, repo_type=RTYPE, operations=ops,
                          commit_message="Rename EC -> EC_v1 (preserve v1 release)"), "rename commit")
        print("    committed.")
    else:
        for f in files[:4]:
            print(f"    EC/{f[3:]}  ->  EC_v1/{f[3:]}  ({'LFS copy' if is_lfs[f] else 're-add'})")
        print(f"    ... ({len(files)} files total)")


def retry(fn, what, tries=5):
    """Retry on transient HF errors (e.g. 504 Gateway Timeout) with backoff."""
    for i in range(tries):
        try:
            return fn()
        except Exception as e:
            if i == tries - 1:
                raise
            wait = 10 * (i + 1)
            print(f"    ! {what} failed ({type(e).__name__}); retry {i+1}/{tries-1} in {wait}s")
            time.sleep(wait)


def upload(api, local: Path, repo_path: str, execute):
    print(f"    {'UPLOAD' if execute else 'would upload'}: {local}  ->  {repo_path}")
    if execute:
        retry(lambda: api.upload_file(path_or_fileobj=str(local), path_in_repo=repo_path,
                                      repo_id=REPO, repo_type=RTYPE), repo_path)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--execute", action="store_true", help="actually push (default: dry run)")
    p.add_argument("--patched-test2", type=Path, default=Path("data_patch/v1_test2/EC"))
    p.add_argument("--v2-root", type=Path, default=Path("v2_build/data/GRIMM_v2"))
    p.add_argument("--card", type=Path, default=Path("data_patch/hf_dataset_card.md"))
    args = p.parse_args()

    api = HfApi()
    mode = "EXECUTE" if args.execute else "DRY RUN (nothing will change)"
    print(f"=== Deploy to {REPO} — {mode} ===\n")

    plan_rename(api, args.execute)

    print("\n[2] overwrite EC_v1 test-2 with decontaminated splits")
    for mod in MODALITIES:
        for s in range(1, 6):
            f = args.patched_test2 / mod / f"split_{s}" / "test2.csv"
            upload(api, f, f"EC_v1/{mod}/split_{s}/test2.csv", args.execute)

    print("\n[3] upload corrected EC_v2 tree (single folder commit)")
    print(f"    {'UPLOAD' if args.execute else 'would upload'}: {args.v2_root}/  ->  EC_v2/")
    if args.execute:
        retry(lambda: api.upload_folder(folder_path=str(args.v2_root), path_in_repo="EC_v2",
                                        repo_id=REPO, repo_type=RTYPE,
                                        commit_message="Add corrected GRIMM-EC v2"), "EC_v2 folder")

    print("\n[4] upload dataset card (README.md)")
    upload(api, args.card, "README.md", args.execute)

    print("\nDone." if args.execute else "\nDry run complete — re-run with --execute to deploy.")


if __name__ == "__main__":
    main()
