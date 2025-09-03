suppressMessages(suppressWarnings({
    library(dada2)
    library(ShortRead)
})

# Get command-line arguments (excluding default R arguments)
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_path  <- args[1]
output_path <- args[2]
path_project_dir <- args[3]

input_path <- "/home/povp/16S_deblur/knihovna11/decontaminated/human_phix"
output_path <- "/home/povp/16S_deblur/knihovna11/dada2/"
path_project_dir <- "/home/povp/16S_deblur/knihovna11"

nproc = 20

# samples
fnFs <- sort(list.files(path, pattern="R1.*\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R1.*\\.fastq\\.gz$", full.names = TRUE))

# initialization
final_seqtab <- NULL

# filtering
sample.names <- sub("_L[0-9]{3}.*", "", basename(files))
filtFs <- file.path(path_to_results, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_to_results, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names
  
# filtering function
out <- filterAndTrim(my_fnFs, filtFs, my_fnRs, filtRs, minLen=200,
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=nproc) 
    
# errors
errF <- learnErrors(filtFs, multithread=nproc,nbases=1e9)
errR <- learnErrors(filtRs, multithread=nproc,nbases=1e9)
    
# dada
dadaFs <- dada(filtFs, err=errF,multithread=nproc)
dadaRs <- dada(filtRs, err=errR, multithread=nproc)
    
# merging pair-end reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, maxMismatch=1, minOverlap=12) 
    
# seqtab
seqtab <- makeSequenceTable(mergers)
    
# chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=nproc, verbose=FALSE)
    
# tracking reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
    
final_seqtab <- seqtab.nochim
final_track <- track

# assigning taxonomy
#taxa <- assignTaxonomy(final_seqtab, paste0(path_to_silva_db,"/silva_nr99_v138.1_train_set.fa.gz"), multithread=10)
#taxa <- addSpecies(taxa, paste0(path_to_silva_db,"silva_species_assignment_v138.1.fa.gz"))

# saving results - track control, taxonomy and abundances
write.csv(final_track,file=file.path(path_to_results,"track_control.csv"))
write.csv(t(final_seqtab), file=file.path(path_to_results,"asv_table.csv")) 

#write.csv(taxa,file=file.path(path_to_results,"taxa_table.csv")) 

write(paste("dada2", packageVersion("dada2")), 
      file = file.path(Sys.getenv("path_project_dir"), "run_info/tools.txt"), 
      append = TRUE)