suppressMessages(suppressWarnings({
    library(dada2)
}))

args <- commandArgs(trailingOnly = TRUE)

# Helper to extract named argument
get_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0) args[idx + 1] else default
}

rget_args_multi <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  start <- idx + 1
  # collect until next flag (starts with --) or end of args
  end <- start
  while (end <= length(args) && !startsWith(args[end], "--")) {
    end <- end + 1
  }
  args[start:(end - 1)]
}

input   <- rget_args_multi(args, "--input")
nproc       <- as.integer(get_arg(args, "--nproc",        "1"))
truncQ      <- as.integer(get_arg(args, "--truncQ",       "2"))
truncLen <- as.integer(get_arg(args, "--truncLen",  "400"))
maxEE    <- as.double(get_arg(args,  "--maxEE",     "2"))

# samples
fnFs <- input

sample.names <- sub("-mergedpairs.fastq.gz", "", basename(fnFs))

# filtering
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))

names(filtFs) <- sample.names
  
# filtering function
out <- filterAndTrim(fnFs, filtFs, truncLen = truncLen, maxN=0, maxEE=maxEE, truncQ=truncQ, rm.phix=TRUE,
                     compress=TRUE, multithread=nproc) 
    
# errors
errF <- learnErrors(filtFs, multithread=nproc,nbases=1e9)
    
# dada
dadaFs <- dada(filtFs, err=errF,multithread=nproc)
# dada() returns a bare dada-class (not a list) when there is only one sample;
# wrap it so downstream list-assuming code works consistently
if (inherits(dadaFs, "dada")) {
    dadaFs <- setNames(list(dadaFs), sample.names)
}

# seqtab
seqtab <- makeSequenceTable(dadaFs)
    
# chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=nproc, verbose=FALSE)

# TO DO - ORIENTING - MERGING


# saving the ASV table
final_seqtab <- as.data.frame(t(seqtab.nochim))
final_seqtab <- cbind(
  SeqID = rownames(final_seqtab),
  final_seqtab
)
write.table(final_seqtab, file="asv_table.tsv",sep="\t",quote=FALSE,row.names=FALSE) 

# tracking reads through the pipeline
getN <- function(x) sum(getUniques(x))
final_track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(final_track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(final_track) <- sample.names 
final_track <- cbind(
  SampleID = rownames(final_track),
  final_track
)

# saving results - track control, taxonomy and abundances
write.table(final_track,file="track_control.tsv",sep="\t",quote=FALSE,row.names=FALSE)

# saving also the FASTA file - for taxonomic assignment purpose
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0(">",colnames(seqtab.nochim))
fasta_lines <- c(rbind(asv_headers, asv_seqs))
writeLines(fasta_lines, "ASV_sequences.fasta")