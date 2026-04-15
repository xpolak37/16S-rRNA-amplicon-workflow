#!/usr/bin/env python3
"""
metastandard_report.py — Self-contained HTML dashboard for MetaStandard results.

Collects all PNG plots and TSV tables produced by the MetaStandard pipeline
processes (per-method plots, cross-method comparison, diversity, agreement)
and bundles them into a single HTML report with embedded base64 images and
styled tables.

Usage:
    python3 metastandard_report.py --input *.png *.tsv --outdir .

File classification is based on filename patterns:
  - *_barplot.png / *_heatmap.png (not compare_/agreement_) -> per-method plots
  - compare_barplot_*.png   -> cross-method barplots
  - compare_braycurtis_*.png -> cross-method heatmaps
  - consensus_table.tsv     -> consensus table
  - diversity_alpha.png     -> alpha diversity plot
  - diversity_alpha.tsv     -> alpha diversity table
  - diversity_pcoa.png      -> PCoA ordination plot
  - agreement_jaccard_*.png -> Jaccard heatmaps
  - agreement_completeness.png -> completeness plot
  - agreement_completeness.tsv -> completeness table
  - agreement_summary.tsv   -> agreement summary
"""

import argparse
import base64
import html
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def img_to_base64(path: Path) -> str:
    """Read a PNG file and return a base64-encoded data URI."""
    data = path.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return f"data:image/png;base64,{b64}"


def tsv_to_html_table(path: Path, max_rows: int = 200, float_fmt: str = ".4f") -> str:
    """Read a TSV and return an HTML table string."""
    df = pd.read_csv(path, sep="\t")
    if len(df) > max_rows:
        df = df.head(max_rows)
        note = f'<p class="meta">Showing first {max_rows} of {len(df)} rows.</p>'
    else:
        note = ""

    # Format floats
    for col in df.select_dtypes(include=["float64", "float32"]).columns:
        df[col] = df[col].map(lambda x: f"{x:{float_fmt}}" if pd.notna(x) else "")

    header = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, row in df.iterrows():
        cells = "".join(f"<td>{html.escape(str(v))}</td>" for v in row)
        rows.append(f"<tr>{cells}</tr>")

    return note + (
        f'<div class="table-wrap"><table>'
        f'<thead><tr>{header}</tr></thead>'
        f'<tbody>{"".join(rows)}</tbody>'
        f'</table></div>'
    )


def img_block(path: Path, caption: str = "") -> str:
    """Return an HTML block with an embedded base64 image."""
    uri = img_to_base64(path)
    cap = f'<p class="caption">{html.escape(caption)}</p>' if caption else ""
    return f'<div class="img-block"><img src="{uri}" alt="{html.escape(caption)}">{cap}</div>'


# ---------------------------------------------------------------------------
# File classification
# ---------------------------------------------------------------------------

def classify_files(paths: list) -> dict:
    """Sort input files into categories based on filename patterns."""
    cats = {
        "per_method_barplots": [],
        "per_method_heatmaps": [],
        "compare_barplots": [],
        "compare_heatmaps": [],
        "consensus_table": None,
        "diversity_alpha_plot": None,
        "diversity_alpha_table": None,
        "diversity_pcoa_plot": None,
        "agreement_jaccard_plots": [],
        "agreement_completeness_plot": None,
        "agreement_completeness_table": None,
        "agreement_summary": None,
    }

    for p in paths:
        p = Path(p)
        name = p.name

        # Agreement files (check before generic barplot/heatmap patterns)
        if name.startswith("agreement_jaccard_") and name.endswith(".png"):
            cats["agreement_jaccard_plots"].append(p)
        elif name == "agreement_completeness.png":
            cats["agreement_completeness_plot"] = p
        elif name == "agreement_completeness.tsv":
            cats["agreement_completeness_table"] = p
        elif name == "agreement_summary.tsv":
            cats["agreement_summary"] = p

        # Cross-method comparison
        elif name.startswith("compare_barplot_") and name.endswith(".png"):
            cats["compare_barplots"].append(p)
        elif name.startswith("compare_braycurtis_") and name.endswith(".png"):
            cats["compare_heatmaps"].append(p)
        elif name == "consensus_table.tsv":
            cats["consensus_table"] = p

        # Diversity
        elif name == "diversity_alpha.png":
            cats["diversity_alpha_plot"] = p
        elif name == "diversity_alpha.tsv":
            cats["diversity_alpha_table"] = p
        elif name == "diversity_pcoa.png":
            cats["diversity_pcoa_plot"] = p

        # Per-method plots (generic barplot/heatmap that aren't compare_ or agreement_)
        elif name.endswith("_barplot.png"):
            cats["per_method_barplots"].append(p)
        elif name.endswith("_heatmap.png"):
            cats["per_method_heatmaps"].append(p)

    # Sort lists for deterministic order
    for key in ["per_method_barplots", "per_method_heatmaps",
                "compare_barplots", "compare_heatmaps",
                "agreement_jaccard_plots"]:
        cats[key].sort(key=lambda x: x.name)

    return cats


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>MetaStandard Report</title>
<style>
  :root {{
    --fg: #1f2933;
    --muted: #52606d;
    --accent: #0b7285;
    --bg: #f7f9fb;
    --card: #ffffff;
    --border: #d9e2ec;
    --nav-bg: #1f2933;
    --nav-fg: #f0f4f8;
    --nav-hover: #0b7285;
    --nav-width: 240px;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    color: var(--fg);
    background: var(--bg);
    margin: 0;
    padding: 0;
    line-height: 1.5;
  }}

  /* Sidebar navigation */
  nav {{
    position: fixed;
    top: 0;
    left: 0;
    width: var(--nav-width);
    height: 100vh;
    background: var(--nav-bg);
    overflow-y: auto;
    padding: 1.5rem 0;
    z-index: 10;
  }}
  nav h2 {{
    color: var(--nav-fg);
    font-size: 1rem;
    padding: 0 1.2rem;
    margin: 0 0 1rem 0;
    letter-spacing: 0.03em;
  }}
  nav a {{
    display: block;
    color: var(--nav-fg);
    text-decoration: none;
    padding: 0.45rem 1.2rem;
    font-size: 0.88rem;
    opacity: 0.8;
    transition: opacity 0.15s, background 0.15s;
  }}
  nav a:hover {{
    opacity: 1;
    background: rgba(255,255,255,0.08);
  }}
  nav .nav-group {{
    margin-top: 1rem;
    padding-top: 0.5rem;
    border-top: 1px solid rgba(255,255,255,0.1);
  }}
  nav .nav-group-title {{
    font-size: 0.72rem;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    color: rgba(255,255,255,0.45);
    padding: 0 1.2rem;
    margin-bottom: 0.3rem;
  }}

  /* Main content */
  main {{
    margin-left: var(--nav-width);
    padding: 2rem 2.5rem;
    max-width: 1200px;
  }}
  h1 {{
    margin-top: 0;
    color: var(--accent);
    font-size: 1.6rem;
  }}
  h2 {{
    color: var(--accent);
    font-size: 1.2rem;
    margin-top: 2rem;
    margin-bottom: 0.5rem;
    padding-bottom: 0.3rem;
    border-bottom: 2px solid var(--border);
  }}
  h3 {{
    color: var(--fg);
    font-size: 1rem;
    margin-top: 1.2rem;
    margin-bottom: 0.3rem;
  }}
  .meta {{
    color: var(--muted);
    font-size: 0.88rem;
    margin-bottom: 1rem;
  }}
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-bottom: 1.25rem;
    box-shadow: 0 1px 2px rgba(0,0,0,0.03);
  }}

  /* Images */
  .img-block {{
    margin: 0.8rem 0;
    text-align: center;
  }}
  .img-block img {{
    max-width: 100%;
    height: auto;
    border: 1px solid var(--border);
    border-radius: 4px;
  }}
  .caption {{
    color: var(--muted);
    font-size: 0.82rem;
    margin-top: 0.3rem;
  }}

  /* Tables */
  .table-wrap {{
    overflow-x: auto;
    margin: 0.8rem 0;
  }}
  table {{
    border-collapse: collapse;
    width: 100%;
    font-size: 0.85rem;
  }}
  th, td {{
    text-align: left;
    padding: 0.4rem 0.6rem;
    border-bottom: 1px solid var(--border);
    white-space: nowrap;
  }}
  th {{
    background: #f0f4f8;
    position: sticky;
    top: 0;
    z-index: 1;
  }}
  tr:hover td {{
    background: #f7fafc;
  }}

  .empty {{
    color: var(--muted);
    font-style: italic;
    padding: 1rem 0;
  }}

  /* Grid for side-by-side plots */
  .plot-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(450px, 1fr));
    gap: 1rem;
  }}
</style>
</head>
<body>

<nav>
  <h2>MetaStandard</h2>
  <a href="#overview">Overview</a>
  <div class="nav-group">
    <div class="nav-group-title">Per-method</div>
    <a href="#per-method-barplots">Barplots</a>
    <a href="#per-method-heatmaps">Heatmaps</a>
  </div>
  <div class="nav-group">
    <div class="nav-group-title">Cross-method</div>
    <a href="#compare-barplots">Comparison barplots</a>
    <a href="#compare-heatmaps">Bray-Curtis heatmaps</a>
    <a href="#consensus">Consensus table</a>
  </div>
  <div class="nav-group">
    <div class="nav-group-title">Diversity</div>
    <a href="#alpha-diversity">Alpha diversity</a>
    <a href="#beta-diversity">PCoA ordination</a>
  </div>
  <div class="nav-group">
    <div class="nav-group-title">Agreement</div>
    <a href="#jaccard">Jaccard similarity</a>
    <a href="#completeness">Classification completeness</a>
    <a href="#agreement-summary">Summary metrics</a>
  </div>
</nav>

<main>
<h1 id="overview">MetaStandard Report</h1>
<p class="meta">Comparative analysis of taxonomic profiles across denoiser and classifier combinations.</p>

{content}

</main>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def section_per_method(cats: dict) -> str:
    """Build HTML for per-method barplots and heatmaps."""
    parts = []

    # Barplots
    parts.append('<h2 id="per-method-barplots">Per-method barplots</h2>')
    if cats["per_method_barplots"]:
        parts.append('<div class="card"><div class="plot-grid">')
        for p in cats["per_method_barplots"]:
            label = p.stem.replace("_barplot", "").replace("_", " ")
            parts.append(img_block(p, label))
        parts.append('</div></div>')
    else:
        parts.append('<p class="empty">No per-method barplots found.</p>')

    # Heatmaps
    parts.append('<h2 id="per-method-heatmaps">Per-method heatmaps</h2>')
    if cats["per_method_heatmaps"]:
        parts.append('<div class="card"><div class="plot-grid">')
        for p in cats["per_method_heatmaps"]:
            label = p.stem.replace("_heatmap", "").replace("_", " ")
            parts.append(img_block(p, label))
        parts.append('</div></div>')
    else:
        parts.append('<p class="empty">No per-method heatmaps found.</p>')

    return "\n".join(parts)


def section_comparison(cats: dict) -> str:
    """Build HTML for cross-method comparison."""
    parts = []

    parts.append('<h2 id="compare-barplots">Cross-method comparison barplots</h2>')
    if cats["compare_barplots"]:
        parts.append('<div class="card">')
        for p in cats["compare_barplots"]:
            sample = p.stem.replace("compare_barplot_", "")
            parts.append(img_block(p, f"Sample: {sample}"))
        parts.append('</div>')
    else:
        parts.append('<p class="empty">No cross-method barplots found.</p>')

    parts.append('<h2 id="compare-heatmaps">Bray-Curtis distance heatmaps</h2>')
    if cats["compare_heatmaps"]:
        parts.append('<div class="card"><div class="plot-grid">')
        for p in cats["compare_heatmaps"]:
            sample = p.stem.replace("compare_braycurtis_", "")
            parts.append(img_block(p, f"Sample: {sample}"))
        parts.append('</div></div>')
    else:
        parts.append('<p class="empty">No Bray-Curtis heatmaps found.</p>')

    parts.append('<h2 id="consensus">Consensus table</h2>')
    if cats["consensus_table"]:
        parts.append('<div class="card">')
        parts.append('<p class="meta">Taxa detected across methods. Higher total detections = more consistent across pipelines.</p>')
        parts.append(tsv_to_html_table(cats["consensus_table"]))
        parts.append('</div>')
    else:
        parts.append('<p class="empty">No consensus table found.</p>')

    return "\n".join(parts)


def section_diversity(cats: dict) -> str:
    """Build HTML for diversity analysis."""
    parts = []

    parts.append('<h2 id="alpha-diversity">Alpha diversity</h2>')
    if cats["diversity_alpha_plot"]:
        parts.append('<div class="card">')
        parts.append(img_block(cats["diversity_alpha_plot"],
                               "Shannon, Simpson, and observed richness across methods"))
        parts.append('</div>')
    else:
        parts.append('<p class="empty">No alpha diversity plot found.</p>')

    if cats["diversity_alpha_table"]:
        parts.append('<div class="card">')
        parts.append('<h3>Alpha diversity values</h3>')
        parts.append(tsv_to_html_table(cats["diversity_alpha_table"]))
        parts.append('</div>')

    parts.append('<h2 id="beta-diversity">PCoA ordination (Bray-Curtis)</h2>')
    if cats["diversity_pcoa_plot"]:
        parts.append('<div class="card">')
        parts.append(img_block(cats["diversity_pcoa_plot"],
                               "Points close together have similar taxonomic profiles"))
        parts.append('</div>')
    else:
        parts.append('<p class="empty">No PCoA plot found.</p>')

    return "\n".join(parts)


def section_agreement(cats: dict) -> str:
    """Build HTML for classifier agreement metrics."""
    parts = []

    parts.append('<h2 id="jaccard">Jaccard similarity</h2>')
    if cats["agreement_jaccard_plots"]:
        parts.append('<div class="card"><div class="plot-grid">')
        for p in cats["agreement_jaccard_plots"]:
            sample = p.stem.replace("agreement_jaccard_", "")
            parts.append(img_block(p, f"Sample: {sample}"))
        parts.append('</div></div>')
    else:
        parts.append('<p class="empty">No Jaccard heatmaps found.</p>')

    parts.append('<h2 id="completeness">Classification completeness</h2>')
    if cats["agreement_completeness_plot"]:
        parts.append('<div class="card">')
        parts.append(img_block(cats["agreement_completeness_plot"],
                               "Fraction of abundance classified at each taxonomic rank"))
        parts.append('</div>')
    else:
        parts.append('<p class="empty">No completeness plot found.</p>')

    if cats["agreement_completeness_table"]:
        parts.append('<div class="card">')
        parts.append('<h3>Completeness by rank</h3>')
        parts.append(tsv_to_html_table(cats["agreement_completeness_table"]))
        parts.append('</div>')

    parts.append('<h2 id="agreement-summary">Agreement summary</h2>')
    if cats["agreement_summary"]:
        parts.append('<div class="card">')
        parts.append('<p class="meta">Jaccard similarity statistics and classification completeness at genus/species level.</p>')
        parts.append(tsv_to_html_table(cats["agreement_summary"]))
        parts.append('</div>')
    else:
        parts.append('<p class="empty">No agreement summary found.</p>')

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build self-contained HTML dashboard for MetaStandard results")
    parser.add_argument("--input", required=True, nargs="+",
                        help="All PNG and TSV files to include in the report")
    parser.add_argument("--outdir", default=".", type=Path,
                        help="Output directory (default: current directory)")
    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cats = classify_files(args.input)

    content_parts = [
        section_per_method(cats),
        section_comparison(cats),
        section_diversity(cats),
        section_agreement(cats),
    ]

    content = "\n\n".join(content_parts)
    html_doc = HTML_TEMPLATE.format(content=content)

    out_path = args.outdir / "metastandard_report.html"
    out_path.write_text(html_doc)
    print(f"MetaStandard report written to: {out_path}")


if __name__ == "__main__":
    main()
