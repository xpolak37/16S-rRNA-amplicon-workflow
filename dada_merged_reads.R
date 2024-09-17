###############################################################################
# change to your paths!!!!
path <- "/storage/brno2/home/xpolak37/ikem/spracovane_data/knihovna9/merged_oriented_sub/"
path_to_results <- "/storage/brno2/home/xpolak37/ikem/spracovane_data/knihovna9/dada_results/"
path_idtaxa <- "/storage/brno2/home/xpolak37/ikem/dada_db/SILVA_SSU_r138_2019.RData"
path_addspecies <- "/storage/brno2/home/xpolak37/diplomka/dada_db/silva_species_assignment_v138.1.fa.gz"
#####################################################################################

# ARGS
args = commandArgs(trailingOnly=TRUE)

if (length(args)>0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

library(dada2); packageVersion("dada2")
library(openxlsx)
library(ShortRead)
library(DECIPHER)
load(path_idtaxa)

dir.create(path_to_results)

# filtering
fnFs <- sort(list.files(path, pattern="_oriented.fq", full.names = TRUE))

final_seqtab <- NULL
sample.names <- gsub("_oriented.fq","",basename(fnFs))
names(fnFs) <- sample.names

filtFs <- file.path(path_to_results, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names  
out <- filterAndTrim(fnFs, filtFs, maxN=0, maxEE=3, truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=10) 
    
errF <- learnErrors(filtFs, multithread=10,nbases=1e9)

#dada
dadaFs <- dada(filtFs,err=errF,multithread=10)
mergers <- dadaFs

#seqtab
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=9, verbose=TRUE)
    
getN <- function(x) sum(getUniques(x))
 if (length(sample.names)==1){
      track <- cbind(out, getN(dadaFs), getN(mergers), rowSums(seqtab.nochim))
    } else{
      track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
    }
colnames(track) <- c("input", "filtered", "denoisedF", "merged", "nonchim")
rownames(track) <- sample.names
    
if (is.null(final_seqtab)) {
    final_seqtab <- seqtab.nochim
    final_track <- track
    rownames(final_seqtab) <- sample.names
} else {
   rownames(seqtab.nochim) <- sample.names
   final_seqtab <- mergeSequenceTables(table1=final_seqtab, table2=seqtab.nochim)
   final_track <- rbind(final_track,track)
  }

 
#  taxa<-assignTaxonomy(final_seqtab, "/storage/brno2/home/xpolak37/diplomka/dada_db/silva_nr99_v138.1_train_set.fa.gz", multithread=9)
#  taxa <- addSpecies(taxa, "/storage/brno2/home/xpolak37/diplomka/dada_db/silva_species_assignment_v138.1.fa.gz")
 
dna <- DNAStringSet(getSequences(final_seqtab))
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=10)
 
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
asv_tax <- t(sapply(tax_info, function(x) {
   m <- match(ranks, x$rank) # returns a position where these arguments match
   taxa <- x$taxon[m] # in a column taxon we have information about taxonomy stored
   taxa[startsWith(taxa, "unclassified_")] <- NA
   taxa
}))

colnames(asv_tax) <- ranks
rownames(asv_tax) <- as.character(dna)	
asv_tax <- asv_tax[,-7]

# add species
taxa <- addSpecies(asv_tax, path_addspecies)

colnames(taxa) <- ranks

write.csv(final_track,file=paste(path_to_results,"/track_control.csv",sep=""))
write.csv(taxa,file=paste(path_to_results,"/taxa.csv",sep=""))
write.csv(t(final_seqtab),file=paste(path_to_results,"/final_seqtab.csv",sep=""))



