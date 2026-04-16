#!/usr/bin/env python3
"""
read_tracking.py — Build a read attrition table and plot from pipeline logs.

Parses per-sample read counts from each pipeline stage and produces:
  1. read_tracking.tsv — samples × stages attrition table
  2. read_tracking.png — stacked bar chart showing read fate per sample

Supported log types (auto-detected by filename pattern):
  - *_cutadapt.json          — Cutadapt JSON reports
  - *_host_removal.log       — Bowtie2 host removal stderr
  - *_phix_removal.log       — Bowtie2 PhiX removal stderr
  - track_control.tsv        — DADA2 paired/single tracking
  - unoise_track_control.tsv — VSEARCH UNOISE3 tracking
  - stats.csv                — Deblur per-sample stats
"""

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Parsers for each log type
# ---------------------------------------------------------------------------

def parse_cutadapt_json(path: Path) -> dict:
    """Parse a cutadapt JSON report. Returns {sample_id: {input, output}}."""
    data = json.loads(path.read_text())
    # Sample ID from filename: {sample}_cutadapt.json
    sample = path.stem.replace("_cutadapt", "")

    read_counts = data.get("read_counts", {})
    input_count = read_counts.get("input", 0)
    output_count = read_counts.get("output", 0)
    # For paired-end, cutadapt reports read pairs
    if data.get("paired", False):
        input_count *= 2
        output_count *= 2

    return {sample: {"raw": input_count, "cutadapt": output_count}}


def parse_bowtie2_log(path: Path, stage: str) -> dict:
    """Parse a bowtie2 stderr log. Returns {sample_id: {stage: reads_passing}}.

    Bowtie2 log format:
      N reads; of these:
      ...
      X (Y%) aligned concordantly exactly 1 time
      ...
      Z (W%) overall alignment rate
    The reads that did NOT align are the ones that pass through (--un-conc-gz).
    """
    text = path.read_text()
    sample = path.stem.replace(f"_{stage}", "")

    total_match = re.search(r"(\d+) reads; of these:", text)
    aligned_0_match = re.search(r"(\d+) \([\d.]+%\) aligned concordantly 0 times", text)

    if total_match and aligned_0_match:
        total = int(total_match.group(1))
        unaligned = int(aligned_0_match.group(1))
        # Bowtie2 reports read pairs; multiply by 2 for individual reads
        return {sample: {stage: unaligned * 2, f"pre_{stage}": total * 2}}

    return {}


def parse_dada2_track(path: Path) -> dict:
    """Parse DADA2 track_control.tsv. Returns {sample_id: {stage: count, ...}}."""
    df = pd.read_csv(path, sep="\t")
    results = {}
    for _, row in df.iterrows():
        sample = str(row.get("SampleID", row.iloc[0]))
        entry = {}
        if "input" in row.index:
            entry["denoiser_input"] = int(row["input"])
        if "filtered" in row.index:
            entry["denoiser_filtered"] = int(row["filtered"])
        if "nonchim" in row.index:
            entry["denoiser_final"] = int(row["nonchim"])
        elif "merged" in row.index:
            entry["denoiser_final"] = int(row["merged"])
        results[sample] = entry
    return results


def parse_unoise_track(path: Path) -> dict:
    """Parse unoise_track_control.tsv. Returns {sample_id: {stage: count, ...}}."""
    df = pd.read_csv(path, sep="\t")
    results = {}
    for _, row in df.iterrows():
        sample = str(row["SampleID"])
        results[sample] = {
            "denoiser_input": int(row["input"]),
            "denoiser_filtered": int(row["filtered"]),
            "denoiser_final": int(row["mapped"]),
        }
    return results


def parse_deblur_stats(path: Path) -> dict:
    """Parse Deblur stats.csv. Returns {sample_id: {stage: count, ...}}."""
    df = pd.read_csv(path, sep=",")
    results = {}
    for _, row in df.iterrows():
        sample = str(row.iloc[0])
        entry = {}
        # Deblur stats columns vary; look for common ones
        if "reads-raw" in row.index:
            entry["denoiser_input"] = int(row["reads-raw"])
        if "reads-hit-artifact" in row.index and "reads-raw" in row.index:
            entry["denoiser_filtered"] = int(row["reads-raw"]) - int(row["reads-hit-artifact"])
        if "reads-deblur" in row.index:
            entry["denoiser_final"] = int(row["reads-deblur"])
        elif "reads-hit-reference" in row.index:
            entry["denoiser_final"] = int(row["reads-hit-reference"])
        results[sample] = entry
    return results


# ---------------------------------------------------------------------------
# File classification and aggregation
# ---------------------------------------------------------------------------

STAGES = ["raw", "cutadapt", "host_removal", "phix_removal",
          "denoiser_input", "denoiser_filtered", "denoiser_final"]

STAGE_LABELS = {
    "raw": "Raw reads",
    "cutadapt": "After trimming",
    "host_removal": "After host removal",
    "phix_removal": "After PhiX removal",
    "denoiser_input": "Denoiser input",
    "denoiser_filtered": "Denoiser filtered",
    "denoiser_final": "Denoiser final",
}


def classify_and_parse(files: list) -> dict:
    """Parse all input files and merge into {sample: {stage: count}}."""
    merged = {}

    for f in files:
        f = Path(f)
        name = f.name
        parsed = {}

        if name.endswith("_cutadapt.json"):
            parsed = parse_cutadapt_json(f)
        elif name.endswith("_host_removal.log"):
            parsed = parse_bowtie2_log(f, "host_removal")
        elif name.endswith("_phix_removal.log"):
            parsed = parse_bowtie2_log(f, "phix_removal")
        elif name == "unoise_track_control.tsv":
            parsed = parse_unoise_track(f)
        elif name == "track_control.tsv":
            parsed = parse_dada2_track(f)
        elif name == "stats.csv":
            parsed = parse_deblur_stats(f)
        else:
            print(f"WARNING: unrecognized file {f}, skipping", file=sys.stderr)
            continue

        # Merge into main dict
        for sample, stages in parsed.items():
            if sample not in merged:
                merged[sample] = {}
            merged[sample].update(stages)

    return merged


# ---------------------------------------------------------------------------
# Output: TSV table
# ---------------------------------------------------------------------------

def build_tracking_table(merged: dict, outdir: Path) -> pd.DataFrame:
    """Build and save a read tracking TSV."""
    rows = []
    for sample in sorted(merged.keys()):
        row = {"Sample": sample}
        for stage in STAGES:
            row[STAGE_LABELS[stage]] = merged[sample].get(stage, "")
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "read_tracking.tsv", sep="\t", index=False)
    return df


# ---------------------------------------------------------------------------
# Output: stacked bar chart showing read fate
# ---------------------------------------------------------------------------

def plot_attrition(df: pd.DataFrame, outdir: Path):
    """Bar chart showing reads retained at each stage per sample.

    For large sample counts (>40), switches to a line plot with samples on
    the x-axis and one line per pipeline stage, which remains readable even
    with hundreds of samples.
    """
    samples = df["Sample"].tolist()
    stage_cols = [c for c in df.columns if c != "Sample"]
    n_samples = len(samples)
    n_stages = len(stage_cols)

    # Convert to numeric, fill missing with NaN
    plot_df = df[stage_cols].apply(pd.to_numeric, errors="coerce")

    colors = plt.cm.viridis(np.linspace(0.2, 0.9, n_stages))

    if n_samples <= 40:
        # Grouped bar chart — works well for moderate sample counts
        fig, ax = plt.subplots(figsize=(max(8, n_samples * 0.6), 6))

        x = np.arange(n_samples)
        width = 0.7 / n_stages
        for i, (col, color) in enumerate(zip(stage_cols, colors)):
            values = plot_df[col].fillna(0).values
            offset = (i - n_stages / 2 + 0.5) * width
            ax.bar(x + offset, values, width, label=col, color=color,
                   edgecolor="white", linewidth=0.3)

        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    else:
        # Line plot — scales to hundreds of samples
        fig, ax = plt.subplots(figsize=(max(10, n_samples * 0.06), 6))

        x = np.arange(n_samples)
        for i, (col, color) in enumerate(zip(stage_cols, colors)):
            values = plot_df[col].fillna(np.nan).values
            ax.plot(x, values, color=color, label=col, linewidth=1.2,
                    marker="." if n_samples <= 100 else None, markersize=3)

        # Show a subset of x-tick labels to avoid overlap
        tick_step = max(1, n_samples // 30)
        ax.set_xticks(x[::tick_step])
        ax.set_xticklabels([samples[i] for i in range(0, n_samples, tick_step)],
                           rotation=45, ha="right", fontsize=7)

    ax.set_ylabel("Read count")
    ax.set_title(f"Read attrition across pipeline stages ({n_samples} samples)")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7, frameon=False)

    plt.tight_layout()
    fig.savefig(outdir / "read_tracking.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build read attrition table from pipeline logs")
    parser.add_argument("--input", required=True, nargs="+",
                        help="Log/stats files from pipeline stages")
    parser.add_argument("--outdir", default=".", type=Path,
                        help="Output directory (default: current directory)")
    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    merged = classify_and_parse(args.input)
    if not merged:
        sys.exit("ERROR: no read counts could be extracted from input files")

    df = build_tracking_table(merged, args.outdir)
    plot_attrition(df, args.outdir)

    print(f"Read tracking complete. {len(merged)} samples, outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
