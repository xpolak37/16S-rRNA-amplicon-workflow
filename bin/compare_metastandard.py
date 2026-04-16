#!/usr/bin/env python3
"""
compare_metastandard.py — Cross-method comparison of MetaStandard16S outputs.

Takes all MetaStandard TSVs (one per denoiser×classifier combination) and
produces three outputs:

1. Grouped barplot   — For each sample, top N taxa with bars grouped by method.
2. Bray-Curtis heatmap — Per-sample symmetric heatmap of method dissimilarity.
3. Consensus table   — TSV counting how many methods detect each taxon above a
                        threshold, highlighting universal vs. singleton taxa.

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
from scipy.spatial.distance import braycurtis


# Rank order (same as plot_metastandard.py)
RANKS = ["d", "p", "c", "o", "f", "g", "s"]

# Maximum number of samples for which individual per-sample plots are produced.
# Beyond this, only the first MAX_PLOT_SAMPLES are plotted and a note is printed.
MAX_PLOT_SAMPLES = 20


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def method_label(filepath: Path) -> str:
    """Extract 'denoiser_classifier' from filename like dada2PE_NaiveBayes_run01_genus.tsv."""
    parts = filepath.stem.split("_")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return filepath.stem


def extract_base_label(taxid: str) -> str:
    """Return the deepest non-Unclassified rank label, e.g. 'g__Bacillus'."""
    tax_part = taxid.split("|")[0]
    parts = [p.strip() for p in tax_part.split(";")]

    rank_map = {}
    for part in parts:
        if "__" in part:
            prefix, value = part.split("__", 1)
            rank_map[prefix.strip()] = value.strip()

    for rank in reversed(RANKS):
        value = rank_map.get(rank, "Unclassified")
        if value and value != "Unclassified":
            return f"{rank}__{value}"

    return "d__Unclassified"


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


def build_long_df(method_data: dict) -> pd.DataFrame:
    """Build a long-format DataFrame: method, sample, taxon_label, abundance."""
    rows = []
    for method, df in method_data.items():
        sample_cols = [c for c in df.columns if c != "TaxID"]
        labels = [extract_base_label(t) for t in df["TaxID"]]
        plot_df = df.copy()
        plot_df["label"] = labels
        # Aggregate duplicate labels within a method (can happen at ASV level)
        agg = plot_df.groupby("label")[sample_cols].sum()
        for sample in sample_cols:
            for taxon, abund in agg[sample].items():
                rows.append({
                    "method": method,
                    "sample": sample,
                    "taxon": taxon,
                    "abundance": abund,
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plot: grouped barplot
# ---------------------------------------------------------------------------

def plot_grouped_barplot(long_df: pd.DataFrame, top_n: int, outdir: Path):
    """One grouped barplot per sample — top N taxa, bars grouped by method."""
    methods = sorted(long_df["method"].unique())
    samples = sorted(long_df["sample"].unique())

    if len(samples) > MAX_PLOT_SAMPLES:
        print(f"Note: {len(samples)} samples detected; plotting first "
              f"{MAX_PLOT_SAMPLES} of {len(samples)} for cross-method barplots.",
              file=sys.stderr)
        samples = samples[:MAX_PLOT_SAMPLES]

    # Global top N taxa by mean abundance across all methods and samples
    mean_abund = long_df.groupby("taxon")["abundance"].mean().sort_values(ascending=False)
    top_taxa = mean_abund.index[:top_n].tolist()

    palette = sns.color_palette("tab20", len(methods))
    method_colors = {m: palette[i] for i, m in enumerate(methods)}

    for sample in samples:
        sample_df = long_df[long_df["sample"] == sample].copy()

        # Filter to top taxa and aggregate "Other"
        plot_data = []
        for method in methods:
            mdf = sample_df[sample_df["method"] == method]
            for taxon in top_taxa:
                row = mdf[mdf["taxon"] == taxon]
                abund = row["abundance"].sum() if len(row) > 0 else 0.0
                plot_data.append({"method": method, "taxon": taxon, "abundance": abund})
            # Other
            other_abund = mdf[~mdf["taxon"].isin(top_taxa)]["abundance"].sum()
            if other_abund > 0:
                plot_data.append({"method": method, "taxon": "Other", "abundance": other_abund})

        pdf = pd.DataFrame(plot_data)
        taxa_order = top_taxa + (["Other"] if pdf[pdf["taxon"] == "Other"]["abundance"].sum() > 0 else [])

        fig, ax = plt.subplots(figsize=(max(8, len(taxa_order) * 1.0), 6))
        n_methods = len(methods)
        n_taxa = len(taxa_order)
        bar_width = 0.8 / n_methods
        x = np.arange(n_taxa)

        for i, method in enumerate(methods):
            mdf = pdf[pdf["method"] == method].set_index("taxon")
            heights = [mdf.loc[t, "abundance"] if t in mdf.index else 0.0 for t in taxa_order]
            offset = (i - n_methods / 2 + 0.5) * bar_width
            ax.bar(x + offset, heights, bar_width, label=method,
                   color=method_colors[method], edgecolor="white", linewidth=0.3)

        ax.set_xticks(x)
        ax.set_xticklabels(taxa_order, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Relative abundance")
        ax.set_title(f"Cross-method comparison — {sample} (top {top_n})")
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)

        plt.tight_layout()
        fig.savefig(outdir / f"compare_barplot_{sample}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Plot: Bray-Curtis method agreement heatmap
# ---------------------------------------------------------------------------

def plot_bray_curtis_heatmap(method_data: dict, outdir: Path):
    """Per-sample symmetric heatmap of Bray-Curtis distances between methods."""
    methods = sorted(method_data.keys())
    if len(methods) < 2:
        print("Skipping Bray-Curtis heatmap: fewer than 2 methods.", file=sys.stderr)
        return

    # Build unified taxon list across all methods
    all_taxa = set()
    for df in method_data.values():
        all_taxa.update(df["TaxID"].tolist())
    all_taxa = sorted(all_taxa)

    # Get sample columns from first method
    first_df = next(iter(method_data.values()))
    sample_cols = [c for c in first_df.columns if c != "TaxID"]

    plot_samples = sample_cols
    if len(sample_cols) > MAX_PLOT_SAMPLES:
        print(f"Note: {len(sample_cols)} samples detected; plotting first "
              f"{MAX_PLOT_SAMPLES} of {len(sample_cols)} for Bray-Curtis heatmaps.",
              file=sys.stderr)
        plot_samples = sample_cols[:MAX_PLOT_SAMPLES]

    for sample in plot_samples:
        # Build abundance vectors (taxa × methods), aligned to all_taxa
        vectors = {}
        for method in methods:
            df = method_data[method]
            series = df.set_index("TaxID").reindex(all_taxa, fill_value=0.0)
            if sample in series.columns:
                vectors[method] = series[sample].values.astype(float)
            else:
                vectors[method] = np.zeros(len(all_taxa))

        # Compute pairwise Bray-Curtis
        n = len(methods)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                vi, vj = vectors[methods[i]], vectors[methods[j]]
                # Guard against zero vectors (braycurtis returns NaN for all-zeros)
                if vi.sum() == 0 or vj.sum() == 0:
                    d = 1.0  # maximally dissimilar if one side has no data
                else:
                    d = braycurtis(vi, vj)
                dist_matrix[i, j] = d
                dist_matrix[j, i] = d

        dist_df = pd.DataFrame(dist_matrix, index=methods, columns=methods)

        fig_size = max(5, n * 0.8)
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        sns.heatmap(dist_df, annot=True, fmt=".3f", cmap="YlOrRd",
                    square=True, linewidths=0.5, ax=ax,
                    cbar_kws={"label": "Bray-Curtis distance"})
        ax.set_title(f"Method agreement — {sample}")

        plt.tight_layout()
        fig.savefig(outdir / f"compare_braycurtis_{sample}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Consensus table
# ---------------------------------------------------------------------------

def build_consensus_table(method_data: dict, threshold: float, outdir: Path):
    """Count how many methods detect each taxon above threshold. Write TSV."""
    methods = sorted(method_data.keys())

    # Get sample columns
    first_df = next(iter(method_data.values()))
    sample_cols = [c for c in first_df.columns if c != "TaxID"]

    # Collect all TaxIDs
    all_taxa = set()
    for df in method_data.values():
        all_taxa.update(df["TaxID"].tolist())
    all_taxa = sorted(all_taxa)

    # Build a presence/absence matrix per method (vectorised, avoids O(n^3) loop)
    # For each method, create a DataFrame indexed by TaxID with boolean columns
    method_presence = {}
    for method in methods:
        df = method_data[method].set_index("TaxID")
        # Reindex to all taxa, fill missing with 0
        df = df.reindex(all_taxa, fill_value=0.0)
        method_presence[method] = df[sample_cols] > threshold

    # Build consensus counts
    results = []
    many_samples = len(sample_cols) > MAX_PLOT_SAMPLES

    for taxid in all_taxa:
        label = extract_base_label(taxid)
        row = {"TaxID": taxid, "label": label}
        total = 0
        for sample in sample_cols:
            detecting = [m for m in methods if method_presence[m].loc[taxid, sample]]
            count = len(detecting)
            total += count
            # Only include per-sample columns when sample count is manageable
            if not many_samples:
                row[f"{sample}_count"] = count
                row[f"{sample}_methods"] = ";".join(detecting) if detecting else "none"
        row["total_detections"] = total
        # With many samples, add aggregate stats instead of per-sample columns
        if many_samples:
            row["mean_methods_per_sample"] = round(total / len(sample_cols), 2)
            row["samples_detected_in"] = sum(
                1 for s in sample_cols
                if any(method_presence[m].loc[taxid, s] for m in methods)
            )
        results.append(row)

    consensus_df = pd.DataFrame(results)
    consensus_df = consensus_df.sort_values("total_detections", ascending=False)

    consensus_df.to_csv(outdir / "consensus_table.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cross-method comparison of MetaStandard16S outputs")
    parser.add_argument("--input", required=True, nargs="+",
                        help="MetaStandard TSV files (one per method)")
    parser.add_argument("--top_n", type=int, default=10,
                        help="Top N taxa for grouped barplot (default: 10)")
    parser.add_argument("--threshold", type=float, default=0.01,
                        help="Abundance threshold for consensus table (default: 0.01)")
    parser.add_argument("--outdir", default=".", type=Path,
                        help="Output directory (default: current directory)")
    return parser.parse_args()


def main():
    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    method_data = load_all_tsvs(args.input)
    if len(method_data) < 2:
        sys.exit("ERROR: need at least 2 MetaStandard TSVs for cross-method comparison")

    long_df = build_long_df(method_data)

    plot_grouped_barplot(long_df, args.top_n, args.outdir)
    plot_bray_curtis_heatmap(method_data, args.outdir)
    build_consensus_table(method_data, args.threshold, args.outdir)

    print(f"Cross-method comparison complete. Outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
