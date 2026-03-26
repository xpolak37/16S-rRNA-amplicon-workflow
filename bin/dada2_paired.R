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

input_R1   <- rget_args_multi(args, "--input_R1")
input_R2   <- rget_args_multi(args, "--input_R2")
nproc       <- as.integer(get_arg(args, "--nproc",        "1"))
truncQ      <- as.integer(get_arg(args, "--truncQ",       "2"))
truncLen_R1 <- as.integer(get_arg(args, "--truncLen_R1",  "250"))
truncLen_R2 <- as.integer(get_arg(args, "--truncLen_R2",  "250"))
maxEE_R1    <- as.double(get_arg(args,  "--maxEE_R1",     "2"))
maxEE_R2    <- as.double(get_arg(args,  "--maxEE_R2",     "2"))
maxMismatch <- as.integer(get_arg(args, "--maxMismatch",  "0"))
minOverlap  <- as.integer(get_arg(args, "--minOverlap",   "12"))

# samples
fnFs <- input_R1
fnRs <- input_R2

sample.names <- sub("_R[12]\\.[a-z_]+\\.fastq\\.gz$", "", basename(fnFs))

# filtering
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filtering function
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(truncLen_R1,truncLen_R2), maxN=0, maxEE=c(maxEE_R1,maxEE_R2), truncQ=truncQ, rm.phix=TRUE,
                     compress=TRUE, multithread=nproc) 
    
# errors
errF <- learnErrors(filtFs, multithread=nproc,nbases=1e9)
errR <- learnErrors(filtRs, multithread=nproc,nbases=1e9)
    
# dada
dadaFs <- dada(filtFs, err=errF,multithread=nproc)
dadaRs <- dada(filtRs, err=errR, multithread=nproc)
    
# merging pair-end reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, maxMismatch=maxMismatch, minOverlap=minOverlap) 
    
# seqtab
seqtab <- makeSequenceTable(mergers)
    
# chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=nproc, verbose=FALSE)

# TO DO - ORIENTING, MERGING

# saving the ASV table
final_seqtab <- as.data.frame(t(seqtab.nochim))
final_seqtab <- cbind(
  SeqID = rownames(final_seqtab),
  final_seqtab
)
write.table(final_seqtab, file="asv_table.tsv",sep="\t",quote=FALSE,row.names=FALSE) 

# tracking reads through the pipeline
getN <- function(x) sum(getUniques(x))
final_track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(final_track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
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