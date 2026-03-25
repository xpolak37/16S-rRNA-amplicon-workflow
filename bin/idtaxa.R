suppressMessages(suppressWarnings({
    library(Biostrings)
    library(DECIPHER)
}))

args           <- commandArgs(trailingOnly = TRUE)
fasta          <- args[1]
classifier     <- args[2]
nproc          <- as.integer(args[3])
denoising_tool <- args[4]

# load classifier
load(classifier)

# load fasta
dna <- readDNAStringSet(fasta)

tax_info <- IdTaxa(test=dna, trainingSet=trainingSet_custom, strand="both", processors=nproc)

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Find the maximum taxonomy depth
max_len <- max(sapply(tax_info, function(x) length(x$taxon[-1])))

# Process and pad each taxonomy vector
asv_tax_df <- as.data.frame(do.call(rbind, lapply(tax_info, function(x) {
    taxa <- x$taxon[-1]
    taxa[startsWith(taxa, "unclassified_")] <- "unassigned"
    taxa[startsWith(taxa, "uncultured")]    <- "unassigned"
    c(taxa, rep("unassigned", max_len - length(taxa)))
})), stringsAsFactors = FALSE)

colnames(asv_tax_df) <- ranks[1:ncol(asv_tax_df)]

# Collapse taxonomy ranks into a single semicolon-separated string
prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")

# Collapse taxonomy ranks into a single semicolon-separated string with prefixes
taxonomy_string <- apply(asv_tax_df, 1, function(x) {
    x <- x[x != "unassigned"]
    paste(paste0(prefixes[1:length(x)], x), collapse = ";")
})
# Get minimum confidence per ASV
confidence_df <- sapply(tax_info, function(x) min(x$confidence))

# Build final dataframe
asv_tax_conf_df <- data.frame(
    SeqID      = as.character(dna),
    Taxonomy   = taxonomy_string,
    Confidence = confidence_df,
    stringsAsFactors = FALSE
)

write.table(asv_tax_conf_df,
            file = paste0(denoising_tool, "_taxa_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)