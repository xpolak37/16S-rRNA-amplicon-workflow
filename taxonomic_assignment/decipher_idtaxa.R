suppressMessages(suppressWarnings({
    library(Biostrings)
    library(DECIPHER)
    library(tibble)
    library(data.table)
}))


# Get command-line arguments (excluding default R arguments)
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_path  <- args[1]
output_path <- args[2]
path_project_dir <- args[3]
decipher_classifier <- args[4]

# track version
cat("\ndecipher_idtaxa.R:\n", 
    file = file.path(path_project_dir,"run_info/tools.txt"), append = TRUE)

# load classifier
load(decipher_classifier)

# load asv table
asv_table <- fread(file.path(input_path,"asv_table.tsv"),
                        sep="\t", header=TRUE,check.names=FALSE)

dna <- DNAStringSet(asv_table$SeqID)
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet_custom, strand="both", processors=20)
 
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Find the maximum taxonomy depth
max_len <- max(sapply(tax_info, function(x) length(x$taxon[-1])))

# Process and pad each taxonomy vector
asv_tax_df <- as.data.frame(do.call(rbind, lapply(tax_info, function(x) {
    taxa <- x$taxon[-1]
    taxa[startsWith(taxa, "unclassified_")] <- "unassigned"
    taxa[startsWith(taxa, "uncultured")] <- "unassigned"
    # pad with NA if shorter
    c(taxa, rep("unassigned", max_len - length(taxa)))
})), stringsAsFactors = FALSE)

# assign ASV IDs as row names
colnames(asv_tax_df) <- ranks
rownames(asv_tax_df) <- as.character(dna)
asv_tax_df <- tibble::rownames_to_column(asv_tax_df, "SeqID")

confidence_df <- sapply(tax_info, function(x) {
    conf <- x$confidence
    conf <- min(conf)
    conf
})

asv_tax_conf_df <- asv_tax_df
asv_tax_conf_df$confidence <- confidence_df

write.table(asv_tax_df,file=file.path(output_path,"taxa_table.csv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(asv_tax_conf_df,file=file.path(output_path,"taxa_table_conf.csv"),sep="\t",row.names=FALSE,quote=FALSE)

# track version
cat(paste("DECIPHER",packageVersion("DECIPHER")), "\n", 
    file = file.path(path_project_dir,"run_info/tools.txt"), append = TRUE)

cat(paste("Classifier", basename(decipher_classifier)), "\n", 
    file = file.path(path_project_dir,"run_info/tools.txt"), append = TRUE)

