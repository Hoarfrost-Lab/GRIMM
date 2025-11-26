"""
Visualize ESM embeddings (PCA/UMAP/t-SNE) and measure distribution shift
between datasets using bootstrapped Sliced Wasserstein distance on:

1. ESM embeddings (representation space)
2. Amino-acid composition (sequence space)
3. Nucleotide composition (DNA space, optional)
4. Correlation analysis between sequence- and ESM-level shifts

Results:
- All plots (.png)
- A JSON file ("stats.json") with mean ± 95% CI for Wasserstein and correlation
"""

# -----------------------
# Imports
# -----------------------
import os
import json
import datetime
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
from sklearn.manifold import TSNE
import umap.umap_ as umap
from collections import Counter
from scipy.stats import wasserstein_distance, pearsonr, spearmanr
import itertools

# Helpers
# -----------------------
def make_output_dir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    print(f"📁 Results directory created: {out_dir}")

def format_esm(a):
    if isinstance(a, dict):
        a = a['mean_representations'][33]
    return a

def seq_to_freq_vector(seq, alphabet):
    seq = seq.upper()
    counts = Counter(seq)
    total = sum(counts.get(ch, 0) for ch in alphabet)
    if total == 0:
        return np.zeros(len(alphabet))
    return np.array([counts.get(ch, 0) / total for ch in alphabet])

def sliced_wasserstein(X, Y, n_projections=100, random_state=None):
    rng = np.random.default_rng(random_state)
    dim = X.shape[1]
    distances = []
    for _ in range(n_projections):
        direction = rng.normal(size=dim)
        direction /= np.linalg.norm(direction)
        X_proj = X @ direction
        Y_proj = Y @ direction
        distances.append(wasserstein_distance(X_proj, Y_proj))
    return np.mean(distances)

# -----------------------
# Bootstrap functions
# -----------------------
def bootstrap_ci(func, X, Y, n_boot=200, alpha=0.05):
    """Bootstrap confidence interval for any distance measure."""
    nX, nY = len(X), len(Y)
    boot_vals = []
    for _ in range(n_boot):
        Xs = X[np.random.randint(0, nX, nX)]
        Ys = Y[np.random.randint(0, nY, nY)]
        boot_vals.append(func(Xs, Ys))
    boot_vals = np.array(boot_vals)
    lower = np.percentile(boot_vals, 100 * (alpha / 2))
    upper = np.percentile(boot_vals, 100 * (1 - alpha / 2))
    return boot_vals.mean(), lower, upper

def bootstrap_corr_ci(x, y, n_boot=200, alpha=0.05):
    """Bootstrap Pearson & Spearman correlation confidence intervals."""
    n = len(x)
    pearson_vals, spearman_vals = [], []
    for _ in range(n_boot):
        idx = np.random.randint(0, n, n)
        xb, yb = x[idx], y[idx]
        r_p, _ = pearsonr(xb, yb)
        r_s, _ = spearmanr(xb, yb)
        pearson_vals.append(r_p)
        spearman_vals.append(r_s)
    def ci(vals):
        vals = np.array(vals)
        return vals.mean(), np.percentile(vals, 100 * (alpha / 2)), np.percentile(vals, 100 * (1 - alpha / 2))
    return {"pearson": ci(pearson_vals), "spearman": ci(spearman_vals)}

# -----------------------
# Data loading
# -----------------------
'''
def load_datasets(dataset_paths, path, emb_out_dir, format_esm_flag=True):
    X_list, labels_list, dataset_names, seq_list, seq_list_nt = [], [], [], [], []
    for idx, curr_path in enumerate(dataset_paths):
        data_name = os.path.splitext(os.path.basename(curr_path))[0]
        data_path = os.path.dirname(curr_path)
        dataset_csv = os.path.join(path, curr_path)
        print(f"📂 Loading dataset: {dataset_csv}")
        df = pd.read_csv(dataset_csv, sep='\t', header=0)

        id_names = df['Entry']
        sequences = df['Sequence'].dropna().tolist()
        seq_list.append(sequences)
        if 'Nucleotide' in df.columns:
            seq_list_nt.append(df['Nucleotide'].dropna().tolist())
        else:
            seq_list_nt.append([])

        dataset_vectors = []
        for filename in tqdm(id_names, desc=f"Loading {data_name}", unit="file"):
            emb_path = os.path.join(path, data_path, emb_out_dir, f"{filename}.pt")
            if not os.path.exists(emb_path):
                continue
            e = torch.load(emb_path)
            if format_esm_flag:
                e = format_esm(e)
            dataset_vectors.append(e.numpy())
        if len(dataset_vectors) == 0:
            print(f"⚠️ Warning: no vectors found for {data_name}")
            continue
        X_list.append(np.vstack(dataset_vectors))
        labels_list.append(np.full(len(dataset_vectors), idx))
        dataset_names.append(data_name.replace('_reformated', ''))
    X_all = np.vstack(X_list)
    labels_all = np.concatenate(labels_list)
    X_all = StandardScaler().fit_transform(X_all)
    print(f"✅ Loaded {len(dataset_names)} datasets, total {X_all.shape[0]} samples.")
    return X_all, labels_all, dataset_names, X_list, seq_list, seq_list_nt
'''

def load_datasets(dataset_paths, path, emb_out_dir, format_esm_flag=True):
    """
    Loads multiple datasets with ESM embeddings, sequences, and EC numbers.

    Returns:
        X_all: np.ndarray (scaled concatenated embeddings)
        labels_all: np.ndarray (dataset index per sample)
        dataset_names: list[str]
        X_list: list[np.ndarray] (per-dataset embeddings)
        seq_list: list[list[str]] (amino acid sequences)
        seq_list_nt: list[list[str]] (nucleotide sequences, if present)
        y_all: np.ndarray (encoded EC numbers across all datasets)
        y_list: list[np.ndarray] (EC numbers per dataset)
    """
    X_list, labels_list, dataset_names, seq_list, seq_list_nt, y_list = [], [], [], [], [], []

    for idx, curr_path in enumerate(dataset_paths):
        data_name = os.path.splitext(os.path.basename(curr_path))[0]
        data_path = os.path.dirname(curr_path)
        dataset_csv = os.path.join(path, curr_path)
        print(f"📂 Loading dataset: {dataset_csv}")
        df = pd.read_csv(dataset_csv, sep='\t', header=0)

        # Required columns
        id_names = df['Entry']
        ec_nums = df['EC number'].fillna("Unknown").to_numpy()

        # Sequences
        seqs = df['Sequence'].dropna().tolist() if 'Sequence' in df.columns else []
        seq_list.append(seqs)

        if 'Nucleotide' in df.columns:
            seq_list_nt.append(df['Nucleotide'].dropna().tolist())
        else:
            seq_list_nt.append([])

        dataset_vectors = []
        for filename in tqdm(id_names, desc=f"Loading {data_name}", unit="file"):
            emb_path = os.path.join(path, data_path, emb_out_dir, f"{filename}.pt")
            if not os.path.exists(emb_path):
                continue
            e = torch.load(emb_path)
            if format_esm_flag:
                e = format_esm(e)
            dataset_vectors.append(e.numpy())

        if len(dataset_vectors) == 0:
            print(f"⚠️ Warning: no vectors found for {data_name}")
            continue

        X = np.vstack(dataset_vectors)
        X_list.append(X)
        labels_list.append(np.full(len(dataset_vectors), idx))
        dataset_names.append(data_name.replace('_reformated', ''))
        y_list.append(ec_nums[:len(X)])  # ensure length alignment

    # Concatenate and standardize
    X_all = np.vstack(X_list)
    X_all = StandardScaler().fit_transform(X_all)
    labels_all = np.concatenate(labels_list)

    # Encode EC numbers
    y_all = np.concatenate(y_list)
    label_encoder = LabelEncoder()
    y_all_encoded = label_encoder.fit_transform(y_all)

    print(f"✅ Loaded {len(dataset_names)} datasets, total {X_all.shape[0]} samples.")
    return X_all, labels_all, dataset_names, X_list, seq_list, seq_list_nt, y_all_encoded, y_list

# -----------------------
# Dimensionality reduction
# -----------------------
def compute_embeddings(X_all):
    print("🔹 Computing PCA...")
    X_pca = PCA(n_components=2).fit_transform(X_all)
    print("🔹 Computing UMAP...")
    X_umap = umap.UMAP(n_components=2, random_state=42, verbose=True).fit_transform(X_all)
    print("🔹 Computing t-SNE...")
    X_tsne = TSNE(n_components=2, random_state=42, init='pca', perplexity=30, verbose=1).fit_transform(X_all)
    return X_pca, X_umap, X_tsne

# -----------------------
# Visualization
# -----------------------
def visualize(X_pca, X_umap, X_tsne, labels_all, dataset_names, out_dir):
    reducers = {"PCA": X_pca, "UMAP": X_umap, "t-SNE": X_tsne}
    markers = ['o', 's', '^', 'P', '*', 'X', 'D', 'v', '<', '>']
    fig, axes = plt.subplots(1, len(reducers), figsize=(18, 6))
    for j, (name, X_red) in enumerate(reducers.items()):
        ax = axes[j]
        marker_cycle = itertools.cycle(markers)
        for idx, dataset_name in enumerate(dataset_names):
            marker = next(marker_cycle)
            mask = labels_all == idx
            ax.scatter(X_red[mask, 0], X_red[mask, 1],
                       label=dataset_name, marker=marker, alpha=0.8)
        ax.set_title(name)
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(loc='best')
    plt.savefig(os.path.join(out_dir, "dimensionality_reduction.png"), bbox_inches='tight', dpi=300)
    plt.close()

'''
def visualize_w_ec(X_pca, X_umap, X_tsne, labels_all, dataset_names, y_all, out_dir, label_names=None):
    """
    Create 2x3 grid:
      Top row  — color/marker by dataset
      Bottom row — color by class label
    """

    reducers = {"PCA": X_pca, "UMAP": X_umap, "t-SNE": X_tsne}
    markers = ['o', 's', '^', 'P', '*', 'X', 'D', 'v', '<', '>']
    cmap = plt.get_cmap("tab10" if len(np.unique(y_all)) <= 10 else "tab20")

    fig, axes = plt.subplots(2, len(reducers), figsize=(18, 12))

    # --- Row 1: colored by dataset ---
    for j, (name, X_red) in enumerate(reducers.items()):
        ax = axes[0, j]
        marker_cycle = itertools.cycle(markers)
        for idx, dataset_name in enumerate(dataset_names):
            marker = next(marker_cycle)
            mask = labels_all == idx
            ax.scatter(X_red[mask, 0], X_red[mask, 1],
                       label=dataset_name, marker=marker, alpha=0.8)
        ax.set_title(f"{name} — by dataset")
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(loc='best', fontsize=8)

    # --- Row 2: colored by class labels ---
    for j, (name, X_red) in enumerate(reducers.items()):
        ax = axes[1, j]
        sc = ax.scatter(X_red[:, 0], X_red[:, 1],
                        c=y_all, cmap=cmap, alpha=0.8)
        ax.set_title(f"{name} — by label")
        ax.set_xticks([]); ax.set_yticks([])
        if label_names is not None:
            # optional colorbar with label names
            cbar = plt.colorbar(sc, ax=ax, ticks=range(len(label_names)))
            cbar.ax.set_yticklabels(label_names)
        else:
            plt.colorbar(sc, ax=ax)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "dimensionality_reduction.png_w_ec.png"), bbox_inches='tight', dpi=300)
    plt.clf()
    plt.close()
'''

def visualize_w_ec(X_pca, X_umap, X_tsne, labels_all, dataset_names, y_all, out_dir, label_names=None, y_raw=None):
    """
    Create 2x3 grid:
      Top row  — color/marker by dataset
      Bottom row — color by EC hierarchy (outer-level grouping)
    """

    reducers = {"PCA": X_pca, "UMAP": X_umap, "t-SNE": X_tsne}
    markers = ['o', 's', '^', 'P', '*', 'X', 'D', 'v', '<', '>']

    # ------------------------------------------------------------
    # Derive hierarchy grouping if raw EC numbers are available
    # ------------------------------------------------------------
    if y_raw is not None:
        # Map each EC number to its first component (outer hierarchy)
        outer_classes = []
        for ec in y_raw:
            if isinstance(ec, str) and ec != "Unknown" and "." in ec:
                outer_classes.append(ec.split(".")[0])
            else:
                outer_classes.append("Unknown")
        outer_classes = np.array(outer_classes)

        # Encode outer hierarchy for coloring
        unique_outer = sorted(np.unique(outer_classes), key=lambda x: (x=="Unknown", x))
        outer_to_idx = {c: i for i, c in enumerate(unique_outer)}
        y_outer = np.array([outer_to_idx[c] for c in outer_classes])

        # Build nice labels for colorbar
        hierarchy_labels = [f"{c}. Class" if c != "Unknown" else "Unknown" for c in unique_outer]
    else:
        y_outer = y_all
        hierarchy_labels = label_names if label_names is not None else None

    # Choose color map
    cmap = plt.get_cmap("tab10" if len(np.unique(y_outer)) <= 10 else "tab20")

    fig, axes = plt.subplots(2, len(reducers), figsize=(18, 12))

    # --- Row 1: colored by dataset ---
    for j, (name, X_red) in enumerate(reducers.items()):
        ax = axes[0, j]
        marker_cycle = itertools.cycle(markers)
        for idx, dataset_name in enumerate(dataset_names):
            marker = next(marker_cycle)
            mask = labels_all == idx
            ax.scatter(X_red[mask, 0], X_red[mask, 1],
                       label=dataset_name, marker=marker, alpha=0.8)
        ax.set_title(f"{name} — by dataset")
        ax.set_xticks([]); ax.set_yticks([])
        ax.legend(loc='best', fontsize=8)

    # --- Row 2: colored by EC hierarchy ---
    for j, (name, X_red) in enumerate(reducers.items()):
        ax = axes[1, j]
        sc = ax.scatter(X_red[:, 0], X_red[:, 1],
                        c=y_outer, cmap=cmap, alpha=0.8)
        ax.set_title(f"{name} — by EC hierarchy")
        ax.set_xticks([]); ax.set_yticks([])
        cbar = plt.colorbar(sc, ax=ax, ticks=range(len(np.unique(y_outer))))
        if hierarchy_labels is not None:
            cbar.ax.set_yticklabels(hierarchy_labels)
        cbar.set_label("EC outer hierarchy", rotation=270, labelpad=15)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "dimensionality_reduction_by_ec_hierarchy.png"), bbox_inches='tight', dpi=300)
    plt.close()

# -----------------------
# Wasserstein computation with bootstrap
# -----------------------
def compute_train_to_test_wasserstein(others_X_list, dataset_names, out_dir, mode="ESM", n_proj=100, n_boot=200):
    train_idx = [i for i, n in enumerate(dataset_names) if "train" in n.lower()]
    train_idx = train_idx[0] if train_idx else 0
    train_name, train_X = dataset_names[train_idx], others_X_list[train_idx]

    results = []
    for i, (X, name) in enumerate(zip(others_X_list, dataset_names)):
        if i == train_idx:
            continue
        mean, lo, hi = bootstrap_ci(lambda a, b: sliced_wasserstein(a, b, n_proj), train_X, X, n_boot=n_boot)
        results.append((train_name, name, mean, lo, hi))
        print(f"{mode} Wasserstein({train_name} → {name}) = {mean:.4f} [{lo:.4f}, {hi:.4f}]")

    df = pd.DataFrame(results, columns=["Train", "Dataset", "Mean", "CI_low", "CI_high"])
    plt.figure(figsize=(6,4))
    plt.bar(df["Dataset"], df["Mean"], yerr=[df["Mean"]-df["CI_low"], df["CI_high"]-df["Mean"]],
            color="teal", alpha=0.8, capsize=5)
    plt.ylabel(f"Sliced Wasserstein ({mode})")
    plt.title(f"{mode}-Level Distribution Shift (Train → Others)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{mode.lower()}_train_to_test_wasserstein.png"), dpi=300)
    plt.close()
    return df

def compute_sequence_level_vectors(seq_list, alphabet):
    seq_freqs = []
    for seqs in seq_list:
        if len(seqs) == 0:
            seq_freqs.append(np.zeros((1, len(alphabet))))
            continue
        freq_vectors = np.array([seq_to_freq_vector(s, alphabet) for s in seqs])
        seq_freqs.append(freq_vectors)
    return seq_freqs

# -----------------------
# Correlation analysis (with bootstrap)
# -----------------------
def compare_wasserstein_correlations(seq_df, esm_df, out_dir, n_boot=200):
    merged = pd.merge(seq_df, esm_df, on="Dataset", suffixes=("_seq", "_esm"))
    if len(merged) == 0:
        print("⚠️ Not enough overlapping datasets for correlation.")
        return None
    x = merged["Mean_seq"].values
    y = merged["Mean_esm"].values
    corrs = bootstrap_corr_ci(x, y, n_boot=n_boot)
    plt.figure(figsize=(5,5))
    sns.regplot(x=x, y=y, color="teal")
    plt.xlabel("Sequence-level Wasserstein")
    plt.ylabel("ESM-level Wasserstein")
    plt.title(f"Correlation of Distribution Shifts\nPearson={corrs['pearson'][0]:.2f} (±95% CI)")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "correlation_seq_vs_esm.png"), dpi=300)
    plt.close()
    return corrs

# -----------------------
# Config
# -----------------------
path = "./data"

#can change this after to be nucleotides and another embedding (e.g. looking glass)
DATASET_PATHS = [
    'Stratified/custom/split_{}/train_reformated.csv',
#    'Stratified/custom/split_{}/train_reformated.csv',
    'Stratified/custom/split_{}/valid_reformated.csv', 
    'Stratified/custom/split_{}/test1_reformated.csv', 
    'Stratified/custom/split_{}/test2_reformated.csv',
    'OM-RGC_v2/OM-RGC-GH-pool_labeled.csv',
    'Stratified/custom/split_{}/GH_valid_reformated.csv',
    'Stratified/custom/split_{}/GH_test_reformated.csv',
    'GH_experimental/GH_experimental_filled.csv',
    'CLEAN/price129.csv', 
    'CLEAN/new392.csv', 
    'CLEAN/halogenase.csv'
]


OUT_DIR='./DS_results/all/split_{}/'
emb_out_dir = "esm_data"
format_esm_flag = True
AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
NT_ALPHABET = list("ACGT")

# -----------------------
# -----------------------
# Main
# -----------------------
if __name__ == "__main__":

    for i in range(5,6):

        out_dir = OUT_DIR.format(i) 
        make_output_dir(out_dir)
        stats = {}
        dataset_paths = [p.format(i) for p in DATASET_PATHS]

        X_all, labels_all, dataset_names, X_list, seq_list, seq_list_nt, y_all, y_list = load_datasets(dataset_paths, path, emb_out_dir, format_esm_flag)
        #X_pca, X_umap, X_tsne = compute_embeddings(X_all)
        #visualize(X_pca, X_umap, X_tsne, labels_all, dataset_names, out_dir)
        #visualize_w_ec(X_pca, X_umap, X_tsne, labels_all, dataset_names, y_all, out_dir) #old without heirarchy grouping
        # y_list is a list of EC numbers per dataset, so flatten it:
        #y_raw = np.concatenate(y_list)
        #visualize_w_ec(X_pca, X_umap, X_tsne, labels_all, dataset_names, y_all, out_dir, y_raw=y_raw)


        esm_df = compute_train_to_test_wasserstein(X_list, dataset_names, out_dir, mode="ESM")
        stats["esm_wasserstein"] = esm_df.set_index("Dataset")[["Mean","CI_low","CI_high"]].to_dict("index")

        seq_freqs = compute_sequence_level_vectors(seq_list, AA_ALPHABET)
        seq_df = compute_train_to_test_wasserstein(seq_freqs, dataset_names, out_dir, mode="AA_Sequence")
        stats["aa_sequence_wasserstein"] = seq_df.set_index("Dataset")[["Mean","CI_low","CI_high"]].to_dict("index")

        #right now ignore nucleotide can add later
        if any(len(s) > 0 for s in seq_list_nt):
            nt_freqs = compute_sequence_level_vectors(seq_list_nt, NT_ALPHABET)
            nt_df = compute_train_to_test_wasserstein(nt_freqs, dataset_names, out_dir, mode="Nucleotide")
            stats["nucleotide_wasserstein"] = nt_df.set_index("Dataset")[["Mean","CI_low","CI_high"]].to_dict("index")

        corr = compare_wasserstein_correlations(
            seq_df.rename(columns={"Mean":"Mean_seq"}),
            esm_df.rename(columns={"Mean":"Mean_esm"}),
            out_dir
        )
        corr=None
        if corr:
            stats["correlation"] = corr

        with open(os.path.join(out_dir, "stats.json"), "w") as f:
            json.dump(stats, f, indent=4)
        
        print(f"✅ Saved all statistics with bootstrapped CIs to {out_dir}/stats.json")

