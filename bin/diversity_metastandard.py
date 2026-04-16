#!/usr/bin/env python3
"""
diversity_metastandard.py — Alpha and beta diversity summary across methods.

Takes all MetaStandard TSVs (one per denoiser×classifier combination) and
produces:

1. Alpha diversity stripplot — Shannon, Simpson, observed richness per sample,
   grouped by method.
2. Alpha diversity table     — TSV with all metrics per method×sample.
3. PCoA ordination plot      — Bray-Curtis PCoA coloured by method, shaped by
   sample.  Shows whether variation is dominated by sample identity or method
   choice.

Input TSV naming convention: {denoiser}_{classifier}_{run_id}_{level}.tsv
Method label is derived as: {denoiser}_{classifier}
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import braycurtis, pdist, squareform


# Rank order (same as other metastandard scripts)
RANKS = ["d", "p", "c", "o", "f", "g", "s"]

# Maximum method×sample combinations for PCoA.  Beyond this, a random subset
# of samples is drawn so the distance matrix stays tractable.
MAX_PCOA_POINTS = 500


# ---------------------------------------------------------------------------
# Helpers (shared pattern with compare_metastandard.py)
# ---------------------------------------------------------------------------

def method_label(filepath: Path) -> str:
    """Extract 'denoiser_classifier' from filename like dada2PE_NaiveBayes_run01_genus.tsv."""
    parts = filepath.stem.split("_")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return filepath.stem


def load_all_tsvs(tsv_paths: list) -> dict:
    """Load all TSVs into {method_label: DataFrame} dict."""
    data = {}
    for p in tsv_paths:
        p = Path(p)
        df = pd.read_csv(p, sep="\t")
        if "TaxID" not in df.columns:
            print(f"WARNING: skipping {p} — no TaxID column", file=sys.stderr)
            continue
        label = method_label(p)
        data[label] = df
    return data


# ---------------------------------------------------------------------------
# Alpha diversity
# ---------------------------------------------------------------------------

def shannon(abundances: np.ndarray) -> float:
    """Shannon entropy (natural log) of a relative abundance vector."""
    p = abundances[abundances > 0]
    return -np.sum(p * np.log(p))


def simpson(abundances: np.ndarray) -> float:
    """Simpson diversity index (1 - D)."""
    return 1.0 - np.sum(abundances ** 2)


def observed_richness(abundances: np.ndarray, threshold: float = 0.0) -> int:
    """Count of taxa with abundance strictly above threshold."""
    return int(np.sum(abundances > threshold))


def compute_alpha_diversity(method_data: dict) -> pd.DataFrame:
    """Compute Shannon, Simpson, observed richness for each method×sample."""
    rows = []
    for method, df in method_data.items():
        sample_cols = [c for c in df.columns if c != "TaxID"]
        for sample in sample_cols:
            abundances = df[sample].values.astype(float)
            # Normalise to relative abundance if not already (sum to 1)
            total = abundances.sum()
            if total > 0:
                rel = abundances / total
            else:
                rel = abundances
            rows.append({
                "method": method,
                "sample": sample,
                "Shannon": shannon(rel),
                "Simpson": simpson(rel),
                "Observed_richness": observed_richness(rel),
            })
    return pd.DataFrame(rows)


def plot_alpha_diversity(alpha_df: pd.DataFrame, outdir: Path):
    """Stripplot + boxplot for each alpha diversity metric, grouped by method."""
    metrics = ["Shannon", "Simpson", "Observed_richness"]
    titles = ["Shannon entropy", "Simpson diversity (1-D)", "Observed richness"]

    n_methods = alpha_df["method"].nunique()
    palette = sns.color_palette("tab10", n_methods)

    fig, axes = plt.subplots(1, 3, figsize=(5 * 3, 5))

    for ax, metric, title in zip(axes, metrics, titles):
        sns.boxplot(data=alpha_df, x="method", y=metric, ax=ax,
                    color="white", fliersize=0, linewidth=0.8)
        sns.stripplot(data=alpha_df, x="method", y=metric, hue="method",
                      palette=palette, size=6, jitter=True, alpha=0.8,
                      legend=False, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel(metric.replace("_", " "))
        ax.tick_params(axis="x", rotation=45)
        for label in ax.get_xticklabels():
            label.set_ha("right")

    plt.tight_layout()
    fig.savefig(outdir / "diversity_alpha.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Beta diversity — PCoA on Bray-Curtis
# ---------------------------------------------------------------------------

def pcoa(dist_matrix: np.ndarray) -> tuple:
    """Classical PCoA (Torgerson) on a square distance matrix.

    Returns (coords, explained_variance_ratio) where coords is (n, 2).
    """
    n = dist_matrix.shape[0]
    # Double-centre the squared distance matrix
    D2 = dist_matrix ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D2 @ H

    # Eigen-decomposition (B is symmetric)
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep only positive eigenvalues
    positive = eigenvalues > 0
    eigenvalues = eigenvalues[positive]
    eigenvectors = eigenvectors[:, positive]

    # Coordinates for first 2 axes
    n_axes = min(2, len(eigenvalues))
    coords = eigenvectors[:, :n_axes] * np.sqrt(eigenvalues[:n_axes])

    total_var = eigenvalues.sum()
    if total_var > 0:
        explained = eigenvalues[:n_axes] / total_var
    else:
        explained = np.zeros(n_axes)

    return coords, explained


def plot_pcoa(method_data: dict, outdir: Path):
    """PCoA ordination on Bray-Curtis distances across all method×sample combinations."""
    methods = sorted(method_data.keys())

    # Build unified taxon list
    all_taxa = set()
    for df in method_data.values():
        all_taxa.update(df["TaxID"].tolist())
    all_taxa = sorted(all_taxa)

    # Get sample columns from first method
    first_df = next(iter(method_data.values()))
    sample_cols = [c for c in first_df.columns if c != "TaxID"]

    # If methods × samples exceeds MAX_PCOA_POINTS, subsample samples
    n_total = len(methods) * len(sample_cols)
    pcoa_samples = sample_cols
    if n_total > MAX_PCOA_POINTS:
        max_per_method = MAX_PCOA_POINTS // len(methods)
        max_per_method = max(2, max_per_method)  # keep at least 2
        rng = np.random.default_rng(42)
        pcoa_samples = sorted(rng.choice(sample_cols, size=min(max_per_method, len(sample_cols)), replace=False))
        print(f"Note: {n_total} method×sample points would be too large for PCoA; "
              f"subsampled to {len(pcoa_samples)} samples ({len(methods) * len(pcoa_samples)} points).",
              file=sys.stderr)

    # Build abundance matrix: rows = method×sample combinations, cols = taxa
    labels_method = []
    labels_sample = []
    matrix_rows = []

    for method in methods:
        df = method_data[method]
        series = df.set_index("TaxID").reindex(all_taxa, fill_value=0.0)
        for sample in pcoa_samples:
            if sample in series.columns:
                vec = series[sample].values.astype(float)
            else:
                vec = np.zeros(len(all_taxa))
            # Normalise
            total = vec.sum()
            if total > 0:
                vec = vec / total
            matrix_rows.append(vec)
            labels_method.append(method)
            labels_sample.append(sample)

    n = len(matrix_rows)
    if n < 3:
        print("Skipping PCoA: fewer than 3 method×sample combinations.", file=sys.stderr)
        return

    # Replace zero vectors with tiny uniform to avoid NaN from braycurtis
    abundance_matrix = np.array(matrix_rows)
    zero_rows = abundance_matrix.sum(axis=1) == 0
    if zero_rows.any():
        abundance_matrix[zero_rows] = 1.0 / abundance_matrix.shape[1]

    # Vectorised pairwise Bray-Curtis (scipy C implementation, much faster)
    condensed = pdist(abundance_matrix, metric="braycurtis")

    # Guard: replace any remaining NaN with 1.0
    condensed = np.nan_to_num(condensed, nan=1.0)

    dist_matrix = squareform(condensed)

    coords, explained = pcoa(dist_matrix)

    if coords.shape[1] < 2:
        print("Skipping PCoA plot: fewer than 2 positive eigenvalues.", file=sys.stderr)
        return

    # Plot
    unique_methods = sorted(set(labels_method))
    unique_samples = sorted(set(labels_sample))

    method_colors = {m: c for m, c in zip(unique_methods, sns.color_palette("tab10", len(unique_methods)))}

    # With many samples, don't use per-sample markers (too many shapes)
    use_sample_markers = len(unique_samples) <= 12
    markers = ["o", "s", "D", "^", "v", "P", "*", "X", "h", "<", ">", "p"]
    if use_sample_markers:
        sample_markers = {s: markers[i % len(markers)] for i, s in enumerate(unique_samples)}
    else:
        sample_markers = {s: "o" for s in unique_samples}

    fig, ax = plt.subplots(figsize=(8, 7))

    for i in range(n):
        ax.scatter(
            coords[i, 0], coords[i, 1],
            color=method_colors[labels_method[i]],
            marker=sample_markers[labels_sample[i]],
            s=60 if len(unique_samples) > 20 else 90,
            edgecolors="black", linewidth=0.3, zorder=3,
            alpha=0.7 if len(unique_samples) > 20 else 1.0,
        )

    # Method legend (colour) — always shown
    method_handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                    markerfacecolor=method_colors[m], markersize=8, label=m)
        for m in unique_methods
    ]
    leg1 = ax.legend(handles=method_handles, title="Method",
                     bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    ax.add_artist(leg1)

    # Sample legend (shape) — only when few enough samples to be meaningful
    if use_sample_markers:
        sample_handles = [
            plt.Line2D([0], [0], marker=sample_markers[s], color="w",
                        markerfacecolor="grey", markersize=8, label=s)
            for s in unique_samples
        ]
        ax.legend(handles=sample_handles, title="Sample",
                  bbox_to_anchor=(1.01, 0.5), loc="center left", fontsize=8, frameon=False)

    pct1 = explained[0] * 100
    pct2 = explained[1] * 100
    ax.set_xlabel(f"PCoA 1 ({pct1:.1f}%)")
    ax.set_ylabel(f"PCoA 2 ({pct2:.1f}%)")
    subtitle = ""
    if len(pcoa_samples) < len(sample_cols):
        subtitle = f" ({len(pcoa_samples)}/{len(sample_cols)} samples shown)"
    ax.set_title(f"Bray-Curtis PCoA — methods × samples{subtitle}")
    ax.axhline(0, color="grey", linewidth=0.3, zorder=1)
    ax.axvline(0, color="grey", linewidth=0.3, zorder=1)

    plt.tight_layout()
    fig.savefig(outdir / "diversity_pcoa.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Alpha and beta diversity summary across MetaStandard methods")
    parser.add_argument("--input", required=True, nargs="+",
                        help="MetaStandard TSV files (one per method)")
    parser.add_argument("--outdir", default=".", type=Path,
                        help="Output directory (default: current directory)")
    return parser.parse_args()


def main():
    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    method_data = load_all_tsvs(args.input)
    if len(method_data) < 2:
        sys.exit("ERROR: need at least 2 MetaStandard TSVs for diversity comparison")

    # Alpha diversity
    alpha_df = compute_alpha_diversity(method_data)
    alpha_df.to_csv(args.outdir / "diversity_alpha.tsv", sep="\t", index=False)
    plot_alpha_diversity(alpha_df, args.outdir)

    # Beta diversity — PCoA
    plot_pcoa(method_data, args.outdir)

    print(f"Diversity summary complete. Outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
