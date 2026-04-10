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
              Supported: domain, phylum, class, order, family, genus, species
              Default: genus

--run_id      A label appended to the output filename to track run parameters.
              Default: run01

OUTPUT
------
A tab-separated file named:  <tool>_<taxa_tool>_<run_id>_<level>.tsv
Written to the current working directory.

Columns:
  SeqID    - Full taxonomy string in semicolon-separated format up to requested level
             e.g. d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus
  <sample> - Raw counts summed across all ASVs assigned to that taxon

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
        help="Taxonomic level (species, genus, family...)"
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


def merge_taxa_genus(asv_table, taxa_table):
    
    # Parse taxonomy column into separate rank columns
    tax_parsed = taxa_table["Taxonomy"].apply(parse_taxonomy)
    tax_df = pd.DataFrame(tax_parsed.tolist(), index=taxa_table["SeqID"])
    
    # Keep only up to genus level
    rank_prefixes = ["d", "p", "c", "o", "f", "g"]
    tax_df = tax_df.reindex(columns=rank_prefixes, fill_value="Unclassified")
    
    # Replace empty strings, whitespace-only, and NaN with Unclassified
    tax_df = tax_df.replace(r'^\s*$', "Unclassified", regex=True)
    tax_df = tax_df.fillna("Unclassified")
    
    # Merge taxonomy with ASV table on SeqID
    asv_indexed = asv_table.set_index("SeqID")
    merged = tax_df.join(asv_indexed, how="outer")
    
    # Get sample columns (everything after taxonomy rank columns)
    sample_cols = [col for col in merged.columns if col not in rank_prefixes]
    
    # Fill taxonomy columns with "Unclassified" (not 0) so that ASVs with
    # no taxonomy match don't produce spurious "d__0;p__0;..." TaxID keys.
    # Fill count columns with 0 for taxonomy entries with no ASV counts.
    merged[rank_prefixes] = merged[rank_prefixes].fillna("Unclassified")
    merged[sample_cols] = merged[sample_cols].fillna(0)
    grouped = merged.groupby(rank_prefixes)[sample_cols].sum().reset_index()
    
    # Reconstruct taxonomy string with prefixes
    grouped["SeqID"] = grouped[rank_prefixes].apply(
    lambda row: ";".join([f"{prefix}__{row[prefix]}" 
                          for prefix in rank_prefixes]),  # <-- removed the if condition
    axis=1
)
    
    # Drop individual rank columns and reorder
    grouped = grouped.drop(columns=rank_prefixes)
    grouped = grouped[["SeqID"] + sample_cols]
    grouped = grouped.rename(columns={"SeqID": "TaxID"})

    
    return grouped


def main():
    args = parse_args()
    
    asv_table  = pd.read_csv(args.asv_table,  sep="\t")
    taxa_table = pd.read_csv(args.taxa_table, sep="\t")
    asv_tool = detect_tool(args.taxa_table)
    result_df = merge_taxa_genus(asv_table, taxa_table)
    
    # Convert counts to relative abundances (columns sum to 1)
    sample_cols = [col for col in result_df.columns if col != "TaxID"]
    result_df[sample_cols] = result_df[sample_cols].div(result_df[sample_cols].sum(axis=0), axis=1)
   
    # saving the final merged table
    outfile = f"{asv_tool}_{args.taxa_tool}_{args.run_id}_{args.level}.tsv"
    result_df.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()

