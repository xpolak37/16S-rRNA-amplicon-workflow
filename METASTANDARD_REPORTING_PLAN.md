# MetaStandard Reporting Improvements — Implementation Plan

## Context

This is a Nextflow DSL2 pipeline for 16S rRNA amplicon benchmarking. It runs up to **16 denoiser×classifier combinations** (4 denoisers: DADA2-PE, DADA2-SE, UNOISE3, Deblur × 4 classifiers: NaiveBayes, BLAST, IDTAXA, AssignTaxonomy). Each combination produces a MetaStandard TSV via `bin/metastandard16S.py` — a unified table of relative abundances (taxa × samples). Currently `bin/plot_metastandard.py` produces a stacked barplot and heatmap **per TSV independently**. There is no cross-method comparison.

### Key files
- `bin/metastandard16S.py` — merges ASV + taxonomy tables into unified TSV (relative abundances)
- `bin/plot_metastandard.py` — produces stacked barplot (top N) + heatmap per TSV
- `modules/MetaStandard16S.nf` — METASTANDARD process definition
- `modules/metastandard_plots.nf` — METASTANDARD_PLOTS process (calls plot_metastandard.py)
- `main.nf` — wiring; all 16 METASTANDARD outputs are collected into `ch_metastandard_plots` via `mix()`, then `METASTANDARD_PLOTS(ch_metastandard_plots)` is called once per TSV
- `nextflow.config` — params including `metastandard_top_n = 10`, `tax_level = "asv"`, `run_id = "run01"`
- `test/mock_test/test_metastandard.tsv` — example MetaStandard output for testing
- Containers: `seaborn:0.13.2` (includes pandas, numpy, scipy, matplotlib) is used for plotting

### TSV naming convention
Output filenames follow: `{denoiser}_{classifier}_{run_id}_{level}.tsv`
e.g., `dada2PE_NaiveBayes_run01_genus.tsv`

---

## Changes to implement (in priority order)

### 1. Cross-method comparison plots (HIGH PRIORITY)

**Goal:** Side-by-side visual comparison of taxonomic profiles across all denoiser×classifier combinations.

**What to build:**
- A new Python script `bin/compare_metastandard.py` that takes **all** MetaStandard TSVs as input (not one at a time)
- **Grouped barplot:** For each sample, show top N taxa with bars grouped by method. X-axis = taxa, color = method. One figure per sample (or a faceted grid for all samples).
- **Method agreement heatmap:** For each sample, compute Bray-Curtis distance between every pair of methods. Display as a symmetric heatmap (methods × methods).
- **Consensus table:** For each taxon, count how many methods detect it above a threshold (e.g., >0.01 relative abundance). Output as TSV + highlight taxa found by all methods vs. only one.

**Nextflow wiring:**
- Currently `ch_metastandard_plots` emits one TSV at a time. You need a new process (e.g., `METASTANDARD_COMPARE`) that takes `ch_metastandard_plots.collect()` — all TSVs at once.
- Add to `modules/metastandard_plots.nf` or create a new module file.
- Call after `METASTANDARD_PLOTS` in `main.nf`.

**Container:** `seaborn:0.13.2` (already available, has scipy for Bray-Curtis).

---

### 2. Alpha/beta diversity summary (MEDIUM PRIORITY)

**Goal:** Show whether method choice or sample identity dominates variation.

**What to build:**
- Extend `bin/compare_metastandard.py` (or a separate script `bin/diversity_metastandard.py`)
- **Alpha diversity:** Compute Shannon, Simpson, observed richness from each MetaStandard TSV. Plot as boxplots/stripplots grouped by method.
- **Beta diversity / ordination:** PCoA or NMDS on Bray-Curtis distances across all samples×methods. Color by method, shape by sample. If samples cluster by sample identity (not method), the methods are consistent.

**Dependencies:** scipy (Bray-Curtis), sklearn or scipy for PCoA. All available in the seaborn container since scipy is a dependency.

---

### 3. Classifier agreement metrics (MEDIUM PRIORITY)

**Goal:** Quantify how much classifiers agree, beyond visual comparison.

**What to build:**
- Per-sample Jaccard similarity of detected taxa (presence/absence above threshold) across all classifiers
- Classification completeness: fraction of reads assigned to each rank (genus, family, etc.) per method
- Summary table output (TSV) + a small heatmap visualization

**Can be part of** `compare_metastandard.py` or a dedicated section in the HTML dashboard (#4).

---

### 4. HTML dashboard (MEDIUM PRIORITY)

**Goal:** Single self-contained HTML report bundling all MetaStandard results.

**What to build:**
- A new script `bin/metastandard_report.py` that reads all TSVs + all generated plots and produces one HTML file
- Sections: per-method barplots/heatmaps (existing), cross-method comparisons (#1), diversity (#2), agreement (#3)
- Sortable/filterable abundance tables (use simple JS or pandas `.to_html()` with classes)
- Embed PNGs as base64 or use inline SVG
- Follow the same pattern as `bin/build_custom_summary.py` (which already produces HTML)

**Nextflow:** A new process `METASTANDARD_REPORT` that takes all plots + TSVs as input and emits `metastandard_report.html`. Runs after all other MetaStandard processes.

---

### 5. Read tracking integration (LOW PRIORITY)

**Goal:** Connect QC metrics to taxonomy in one view — how many reads survived each step per sample.

**What to build:**
- Collect read counts from key stages: raw (FastQC), post-trim (Cutadapt logs), post-host-removal, post-denoising (ASV table column sums)
- Produce an "attrition table" (samples × stages) and a stacked bar or Sankey-style plot
- Include in the HTML dashboard

**Challenge:** Requires wiring outputs from multiple earlier processes into a new collection process. More invasive to `main.nf`.

---

### 6. Differential abundance preview (LOW PRIORITY)

**Goal:** If sample metadata includes group labels, show which taxa differ between groups.

**What to build:**
- Optional: only runs if a grouping column is present in the samplesheet
- Basic Wilcoxon/Kruskal-Wallis test per taxon, volcano plot or MA plot
- Clearly labeled as exploratory, not a formal statistical analysis

**Deferred:** Requires samplesheet metadata support that may not exist yet.

---

## Implementation notes

- All new Python scripts go in `bin/` and should follow the existing pattern (argparse, `if __name__ == "__main__"`)
- Use the `seaborn:0.13.2` container for all plotting processes (already pulled by the pipeline)
- Keep the existing per-TSV plots (`plot_metastandard.py`) — they're still useful individually
- The cross-method comparison (#1) is the most impactful change and should be done first
- Test with `test/mock_test/test_metastandard.tsv` — you can create a few copies with slightly different values to simulate multiple methods
- Publish all outputs to `${params.outdir}/metastandard_report/`
