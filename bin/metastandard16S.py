#!/usr/bin/env python3

"""
MetaStandard 16S - ASV Taxonomic Profile Unifier
=============================================
Merges and standardises ASV-based taxonomic abundance profiles from multiple
denoising tools (DADA2, UNOISE, Deblur) into a single, unified TSV table.

WHAT IT DOES
------------
1. Auto-detects the denoising tool used based on input filenames.
2. Parses ASV and taxonomy tables, linking sequences to taxonomic annotations.
3. Aggregates ASV counts to the requested taxonomic level (default: genus).
   Use --level asv to skip aggregation and keep individual ASV rows.
4. Standardises taxonomy strings using rank prefixes (d__, p__, c__, o__, f__, g__).
5. Merges all samples into one wide-format table (taxa x samples).
6. Writes the result to a TSV file.

INPUTS
------
--asv_table   One or more ASV count tables in TSV format.
              Rows are ASV sequences (SeqID), columns are samples.

--taxa_table  One or more taxonomy tables in TSV format.
              Must contain columns: SeqID, Taxonomy, Confidence
              Taxonomy strings must follow the format:
                d__Bacteria;p__Firmicutes;c__Bacilli;...;g__Lactobacillus

--taxa_tool   Name(s) of the tool(s) used to generate the taxonomy
              (e.g. dada2, unoise, deblur).
              Tool is also auto-detected from filename if not specified.

--level       Taxonomic level to aggregate counts to.
              Supported: domain, phylum, class, order, family, genus, species, asv
              Default: genus
              Use 'asv' to keep individual ASV sequences without aggregation.

--run_id      A label appended to the output filename to track run parameters.
              Default: run01

OUTPUT
------
A tab-separated file named:  <tool>_<taxa_tool>_<run_id>_<level>.tsv
Written to the current working directory.

Columns (taxonomic levels):
  TaxID    - Full taxonomy string in semicolon-separated format up to requested level
             e.g. d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus
  <sample> - Relative abundance summed across all ASVs assigned to that taxon

Columns (asv level):
  TaxID    - Full taxonomy string up to species + "|" + ASV sequence
             e.g. d__Bacteria;...;g__Lactobacillus;s__crispatus|ACGT...
  <sample> - Relative abundance of the individual ASV

Unclassified taxa are labelled as "Unclassified" at each rank.

USAGE EXAMPLES
--------------
# Merge DADA2 ASV and taxonomy tables at genus level:
python metastandard16S.py \\
    --asv_table  asv_table.tsv \\
    --taxa_table dada2_taxa_table.tsv \\
    --taxa_tool  qiime_blast \\
    --level      genus \\
    --run_id     run01

# Keep individual ASV rows with full taxonomy annotation:
python metastandard16S.py \\
    --asv_table  asv_table.tsv \\
    --taxa_table dada2_taxa_table.tsv \\
    --taxa_tool  qiime_blast \\
    --level      asv \\
    --run_id     run01


DEPENDENCIES
------------
  pandas

LIMITATIONS
-----------
- ASV and taxonomy tables must be provided in matching order.
- SeqID column must contain the ASV sequences, not arbitrary IDs.
- Taxonomy strings must use the double-underscore prefix format (d__, p__, etc.).
"""

import argparse
import pandas as pd

RANK_PREFIXES = ["d", "p", "c", "o", "f", "g", "s"]

LEVEL_TO_PREFIX = {
    "domain":  "d",
    "phylum":  "p",
    "class":   "c",
    "order":   "o",
    "family":  "f",
    "genus":   "g",
    "species": "s",
}

def parse_args():
    parser = argparse.ArgumentParser(description="MetaStandard: unify taxonomic profiles")

    parser.add_argument(
        "--asv_table",
        required=True,
        help="Input ASV table"
    )

    parser.add_argument(
        "--taxa_table",
        required=True,
        help="Input taxa table"
    )

    parser.add_argument(
        "--taxa_tool",
        required=True,
        help="Tool used to generate the taxonomy"
    )

    parser.add_argument(
        "--level",
        default="genus",
        help="Taxonomic level (domain, phylum, class, order, family, genus, species, asv)"
    )

    parser.add_argument(
        "--run_id",
        default="run01",
        help="ID of the run in order to recognize the parameters used"
    )

    return parser.parse_args()


def detect_tool(f):

    if "dada2_paired" in f.lower():
        return "dada2PE"

    if "dada2_single" in f.lower():
        return "dada2SE"

    if "unoise" in f.lower():
        return "unoise"

    if "deblur" in f.lower():
        return "deblur"

    return "unknown"


def parse_taxonomy(taxonomy_string):
    """Parse d__Bacteria;p__...;g__Genus style string into a dictionary"""
    ranks = {}
    if pd.isna(taxonomy_string):
        return ranks
    for part in taxonomy_string.split(";"):
        part = part.strip()
        if "__" in part:
            prefix, value = part.split("__", 1)
            ranks[prefix.strip()] = value.strip()
    return ranks


def build_tax_df(taxa_table):
    """Parse taxonomy column into a DataFrame of rank columns (d..s)."""
    tax_parsed = taxa_table["Taxonomy"].apply(parse_taxonomy)
    tax_df = pd.DataFrame(tax_parsed.tolist(), index=taxa_table["SeqID"])
    tax_df = tax_df.reindex(columns=RANK_PREFIXES, fill_value="Unclassified")
    tax_df = tax_df.replace(r'^\s*$', "Unclassified", regex=True)
    tax_df = tax_df.fillna("Unclassified")
    return tax_df


def aggregate_to_level(asv_table, taxa_table, level):
    """Aggregate ASV counts to the requested taxonomic level."""
    depth_prefix = LEVEL_TO_PREFIX[level]
    rank_cols = RANK_PREFIXES[: RANK_PREFIXES.index(depth_prefix) + 1]

    tax_df = build_tax_df(taxa_table)
    tax_df = tax_df.reindex(columns=rank_cols, fill_value="Unclassified")

    asv_indexed = asv_table.set_index("SeqID")
    merged = tax_df.join(asv_indexed, how="outer")

    sample_cols = [col for col in merged.columns if col not in rank_cols]

    merged[rank_cols] = merged[rank_cols].fillna("Unclassified")
    merged[sample_cols] = merged[sample_cols].fillna(0)

    grouped = merged.groupby(rank_cols)[sample_cols].sum().reset_index()

    grouped["TaxID"] = grouped[rank_cols].apply(
        lambda row: ";".join([f"{p}__{row[p]}" for p in rank_cols]),
        axis=1
    )

    grouped = grouped.drop(columns=rank_cols)
    grouped = grouped[["TaxID"] + sample_cols]

    return grouped


def aggregate_to_asv(asv_table, taxa_table):
    """Keep individual ASV rows; TaxID = full taxonomy string | ASV sequence."""
    tax_df = build_tax_df(taxa_table)

    asv_indexed = asv_table.set_index("SeqID")
    merged = tax_df.join(asv_indexed, how="outer")

    sample_cols = [col for col in merged.columns if col not in RANK_PREFIXES]

    merged[RANK_PREFIXES] = merged[RANK_PREFIXES].fillna("Unclassified")
    merged[sample_cols] = merged[sample_cols].fillna(0)

    tax_string = merged[RANK_PREFIXES].apply(
        lambda row: ";".join([f"{p}__{row[p]}" for p in RANK_PREFIXES]),
        axis=1
    )
    merged["TaxID"] = tax_string + "|" + merged.index

    result = merged[["TaxID"] + sample_cols].reset_index(drop=True)

    return result


def main():
    args = parse_args()

    level = args.level.lower()
    if level not in LEVEL_TO_PREFIX and level != "asv":
        raise ValueError(
            f"Unknown level '{args.level}'. "
            f"Choose from: {', '.join(list(LEVEL_TO_PREFIX.keys()) + ['asv'])}"
        )

    asv_table  = pd.read_csv(args.asv_table,  sep="\t")
    taxa_table = pd.read_csv(args.taxa_table, sep="\t")
    asv_tool   = detect_tool(args.taxa_table)

    if level == "asv":
        result_df = aggregate_to_asv(asv_table, taxa_table)
    else:
        result_df = aggregate_to_level(asv_table, taxa_table, level)

    # Convert counts to relative abundances (columns sum to 1)
    sample_cols = [col for col in result_df.columns if col != "TaxID"]
    result_df[sample_cols] = result_df[sample_cols].div(
        result_df[sample_cols].sum(axis=0), axis=1
    )

    outfile = f"{asv_tool}_{args.taxa_tool}_{args.run_id}_{level}.tsv"
    result_df.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
