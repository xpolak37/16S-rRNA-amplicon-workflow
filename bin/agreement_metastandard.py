#!/usr/bin/env python3
"""
agreement_metastandard.py — Classifier agreement metrics across methods.

Takes all MetaStandard TSVs (one per denoiser×classifier combination) and
produces:

1. Jaccard similarity heatmap — Per-sample symmetric heatmap showing
   presence/absence agreement (Jaccard index) between every pair of methods.
2. Classification completeness plot — Fraction of total abundance assigned
   at each taxonomic rank, per method.
3. Agreement summary table — TSV with Jaccard values and completeness metrics.

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


# Rank order (same as other metastandard scripts)
RANKS = ["d", "p", "c", "o", "f", "g", "s"]
RANK_NAMES = {
    "d": "Domain", "p": "Phylum", "c": "Class",
    "o": "Order", "f": "Family", "g": "Genus", "s": "Species",
}


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


def deepest_rank(taxid: str) -> str:
    """Return the deepest classified rank prefix from a TaxID string.

    E.g. 'd__Bacteria;p__Firmicutes;c__Bacilli;o__;f__;g__;s__'
    returns 'c' because class is the deepest non-empty rank.
    A completely unclassified entry returns '' (empty string).
    """
    tax_part = taxid.split("|")[0]
    parts = [p.strip() for p in tax_part.split(";")]
    found = ""
    for part in parts:
        if "__" in part:
            prefix, value = part.split("__", 1)
            prefix = prefix.strip()
            value = value.strip()
            if value and value != "Unclassified":
                found = prefix
    return found


# ---------------------------------------------------------------------------
# Jaccard similarity
# ---------------------------------------------------------------------------

def compute_jaccard_per_sample(method_data: dict, threshold: float) -> dict:
    """Compute pairwise Jaccard similarity of detected taxa for each sample.

    Returns {sample: DataFrame (methods × methods)}.
    """
    methods = sorted(method_data.keys())
    first_df = next(iter(method_data.values()))
    sample_cols = [c for c in first_df.columns if c != "TaxID"]

    results = {}
    for sample in sample_cols:
        # Build detected-taxa sets per method
        detected = {}
        for method in methods:
            df = method_data[method]
            mask = df[sample].values.astype(float) > threshold if sample in df.columns else np.zeros(len(df), dtype=bool)
            detected[method] = set(df.loc[mask, "TaxID"].tolist())

        # Pairwise Jaccard
        n = len(methods)
        matrix = np.ones((n, n))  # diagonal = 1.0
        for i in range(n):
            for j in range(i + 1, n):
                a, b = detected[methods[i]], detected[methods[j]]
                union = len(a | b)
                if union == 0:
                    j_val = 1.0  # both empty → perfect agreement
                else:
                    j_val = len(a & b) / union
                matrix[i, j] = j_val
                matrix[j, i] = j_val

        results[sample] = pd.DataFrame(matrix, index=methods, columns=methods)

    return results


def plot_jaccard_heatmaps(jaccard_data: dict, outdir: Path):
    """One heatmap per sample showing Jaccard similarity between methods."""
    for sample, jdf in jaccard_data.items():
        n = len(jdf)
        fig_size = max(5, n * 0.8)
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        sns.heatmap(jdf, annot=True, fmt=".3f", cmap="YlGnBu",
                    square=True, linewidths=0.5, ax=ax,
                    vmin=0, vmax=1,
                    cbar_kws={"label": "Jaccard similarity"})
        ax.set_title(f"Classifier agreement (Jaccard) — {sample}")

        plt.tight_layout()
        fig.savefig(outdir / f"agreement_jaccard_{sample}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Classification completeness
# ---------------------------------------------------------------------------

def compute_completeness(method_data: dict) -> pd.DataFrame:
    """For each method, compute fraction of total abundance classified at each rank.

    Returns a DataFrame with columns: method, rank, rank_name, fraction.
    """
    rows = []
    for method, df in method_data.items():
        sample_cols = [c for c in df.columns if c != "TaxID"]
        # Total abundance across all samples
        total_abund = df[sample_cols].sum().sum()
        if total_abund == 0:
            continue

        # For each row, determine deepest rank and accumulate abundance
        rank_abund = {r: 0.0 for r in RANKS}
        for _, row in df.iterrows():
            dr = deepest_rank(row["TaxID"])
            row_abund = row[sample_cols].sum()
            if dr in rank_abund:
                rank_abund[dr] += row_abund

        # Cumulative: classified "at least to rank X" = sum of that rank and deeper
        for rank in RANKS:
            # Fraction classified at exactly this depth
            frac = rank_abund[rank] / total_abund if total_abund > 0 else 0.0
            rows.append({
                "method": method,
                "rank": rank,
                "rank_name": RANK_NAMES[rank],
                "fraction": frac,
            })

    return pd.DataFrame(rows)


def compute_cumulative_completeness(completeness_df: pd.DataFrame) -> pd.DataFrame:
    """Convert exact-rank fractions to cumulative: classified at least to rank X.

    E.g. fraction classified to at least genus = genus + species.
    """
    rows = []
    for method in completeness_df["method"].unique():
        mdf = completeness_df[completeness_df["method"] == method]
        rank_fracs = dict(zip(mdf["rank"], mdf["fraction"]))
        for i, rank in enumerate(RANKS):
            # Sum from this rank to the deepest
            cum_frac = sum(rank_fracs.get(r, 0.0) for r in RANKS[i:])
            rows.append({
                "method": method,
                "rank": rank,
                "rank_name": RANK_NAMES[rank],
                "cumulative_fraction": cum_frac,
            })
    return pd.DataFrame(rows)


def plot_completeness(completeness_df: pd.DataFrame, outdir: Path):
    """Grouped bar chart: classification completeness at each rank, per method."""
    cum_df = compute_cumulative_completeness(completeness_df)

    methods = sorted(cum_df["method"].unique())
    n_methods = len(methods)
    palette = sns.color_palette("tab10", n_methods)
    method_colors = {m: palette[i] for i, m in enumerate(methods)}

    rank_labels = [RANK_NAMES[r] for r in RANKS]
    x = np.arange(len(RANKS))
    bar_width = 0.8 / n_methods

    fig, ax = plt.subplots(figsize=(max(8, len(RANKS) * 1.2), 5))

    for i, method in enumerate(methods):
        mdf = cum_df[cum_df["method"] == method].set_index("rank")
        heights = [mdf.loc[r, "cumulative_fraction"] if r in mdf.index else 0.0 for r in RANKS]
        offset = (i - n_methods / 2 + 0.5) * bar_width
        ax.bar(x + offset, heights, bar_width, label=method,
               color=method_colors[method], edgecolor="white", linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(rank_labels, fontsize=9)
    ax.set_ylabel("Fraction of abundance classified")
    ax.set_title("Classification completeness by rank")
    ax.set_ylim(0, 1.05)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    ax.axhline(1.0, color="grey", linewidth=0.3, linestyle="--")

    plt.tight_layout()
    fig.savefig(outdir / "agreement_completeness.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def build_agreement_table(jaccard_data: dict, completeness_df: pd.DataFrame,
                          outdir: Path):
    """Write a combined summary TSV with Jaccard stats and completeness."""
    rows = []

    # Jaccard summary: mean pairwise Jaccard per sample
    for sample, jdf in jaccard_data.items():
        n = len(jdf)
        # Extract upper triangle values (excluding diagonal)
        upper = []
        for i in range(n):
            for j in range(i + 1, n):
                upper.append(jdf.iloc[i, j])
        rows.append({
            "metric": "jaccard_mean",
            "scope": sample,
            "value": np.mean(upper) if upper else np.nan,
        })
        rows.append({
            "metric": "jaccard_min",
            "scope": sample,
            "value": np.min(upper) if upper else np.nan,
        })
        rows.append({
            "metric": "jaccard_max",
            "scope": sample,
            "value": np.max(upper) if upper else np.nan,
        })

    # Completeness: cumulative fraction at genus level per method
    cum_df = compute_cumulative_completeness(completeness_df)
    genus_df = cum_df[cum_df["rank"] == "g"]
    for _, row in genus_df.iterrows():
        rows.append({
            "metric": "completeness_genus",
            "scope": row["method"],
            "value": row["cumulative_fraction"],
        })

    # Completeness: cumulative fraction at species level per method
    species_df = cum_df[cum_df["rank"] == "s"]
    for _, row in species_df.iterrows():
        rows.append({
            "metric": "completeness_species",
            "scope": row["method"],
            "value": row["cumulative_fraction"],
        })

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(outdir / "agreement_summary.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Classifier agreement metrics across MetaStandard methods")
    parser.add_argument("--input", required=True, nargs="+",
                        help="MetaStandard TSV files (one per method)")
    parser.add_argument("--threshold", type=float, default=0.01,
                        help="Abundance threshold for Jaccard presence/absence (default: 0.01)")
    parser.add_argument("--outdir", default=".", type=Path,
                        help="Output directory (default: current directory)")
    return parser.parse_args()


def main():
    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    method_data = load_all_tsvs(args.input)
    if len(method_data) < 2:
        sys.exit("ERROR: need at least 2 MetaStandard TSVs for agreement metrics")

    # Jaccard similarity
    jaccard_data = compute_jaccard_per_sample(method_data, args.threshold)
    plot_jaccard_heatmaps(jaccard_data, args.outdir)

    # Classification completeness
    completeness_df = compute_completeness(method_data)
    completeness_df.to_csv(args.outdir / "agreement_completeness.tsv", sep="\t", index=False)
    plot_completeness(completeness_df, args.outdir)

    # Combined summary table
    build_agreement_table(jaccard_data, completeness_df, args.outdir)

    print(f"Agreement metrics complete. Outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
