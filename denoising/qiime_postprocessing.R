suppressMessages(suppressWarnings({
    library(Biostrings)
}))

# Get command-line arguments (excluding default R arguments)
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_path  <- args[1]
output_path <- args[2]
path_project_dir <- args[3]

path_feature_table <- file.path(input_path,"feature-table.tsv")
system(paste0("sed -i '2s/^#OTU ID/OTU_ID/' ", shQuote(path_feature_table)))

dna_seqs <- readDNAStringSet(file.path(input_path,"dna-sequences.fasta"))
asv_table <- read.table(file.path(input_path,"feature-table.tsv"), 
                        skip=1, header=TRUE, check.names = FALSE)

dna_seqs <- dna_seqs[asv_table$OTU_ID]

asv_table$OTU_ID <- as.character(dna_seqs)

colnames(asv_table)[1] <- "SeqID"

write.table(asv_table,file.path(output_path,"asv_table.tsv"),row.names=FALSE,quote=FALSE,sep="\t")