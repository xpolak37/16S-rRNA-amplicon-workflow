suppressMessages(suppressWarnings({
library("ShortRead")
library("stringr")
library("crayon")
}))

# Get command-line arguments (excluding default R arguments)
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
reads_path  <- args[1]
primers_path <- args[2]
output_path <- args[3]

#' @title cross_contamination: Reveal cross contamination in the experiment
#'
#' @description This function examines raw sequencing data, sample by sample, counting
#' the presence of spacers in the reads. The primers sequences and the
#' reference have to be provided in the additional files. If most counts
#' align with the expected reference spacer, comprising over a certain
#' percentage of the reads (input threshold), the function does not report
#' cross-contamination. If it depathtects a different spacer, indicating that
#' the expected spacer is not dominant or is present in less than a given
#' threshold, it reports a possible contamination.
#'
#' @param reads_path Path to directory where FASTQ files are located
#' @param primers_path path to primer sequences (FASTA), if NULL, 
#' the default primers are loaded located in primers.fasta
#' @param forward_pattern Pattern with which forward reads are recognized, 
#' default "\_L00\\d_R1_00\\d.fastq"
#' @param reverse_pattern Pattern with which reverse reads are recognized,
#' default "\_L00\\d_R2_00\\d.fastq"
#' @param mode  Mode in which cross-contamination status should be generated, 
#' "forward" for forward reads only, "reverse" for reverse reads only, 
#' "both" for both forward and reverse read, default "reverse"
#' @param threshold Proportion of reference counts to be reported as non-contaminated, default 0.9
#' @return cross_contamination_status - a data frame with spacer counts and
#' contamination status for each sample
#' @examples
#' # Example usage:
#' status <- cross_contamination(reads_path="data/run/", primers_path="path/to/primers.fasta")
cross_contamination <- function(reads_path, primers_path=NULL, 
                                forward_pattern="_L00\\d_R1_00\\d.fastq",
                                reverse_pattern="_L00\\d_R2_00\\d.fastq",
                                mode="both", threshold=0.9){
  if (is.null(primers_path)){
    message(blue("Analysing default primers"))
    primers <- readDNAStringSet("../primers.fasta")
  } else primers <- readDNAStringSet(primers_path)  
  
  primersF_names <- names(primers)[grepl("F\\d",names(primers))]
  primersR_names <- names(primers)[grepl("R\\d",names(primers))]
  
  method="row"
  # get sequences
  primersF <- as.character(primers[primersF_names])
  primersR <- as.character(primers[primersR_names])
  
  # prepare primers with IUPAC coding
  iupac_codes <- c("N","M","R","W","S","Y","K", "V","H","D","B")
  iupac_expressions <- c("[ACGT]","[AC]","[AG]","[AT]","[CG]","[CT]","[GT]", "[ACG]","[ACT]","[AGT]","[CGT]")
  
  for (i in 1:length(iupac_codes)) {
    primersF <- str_replace_all(primersF,iupac_codes[i], iupac_expressions[i])
    primersR <- str_replace_all(primersR,iupac_codes[i], iupac_expressions[i])
  }
  
  # get the samples and reads in given path
  if (mode=="forward" | mode=="both") forward_reads <- sort(list.files(paste(reads_path), pattern=forward_pattern, full.names = TRUE))
  if (mode=="reverse" | mode=="both") reverse_reads <- sort(list.files(paste(reads_path), pattern=reverse_pattern, full.names = TRUE))
  # analysis samples one by one in for loop
  if (mode=="forward" | mode=="both") {
    for (num_sample in 1:length(forward_reads)){
      # fastq
      fq <- readFastq(forward_reads[num_sample])
      
      # sequences
      reads <- sread(fq)
      
      # sample
      sample <- c(substring(basename(forward_reads[num_sample]), 1, regexpr(forward_pattern, basename(forward_reads[num_sample]))-1))
      
      # look at each primer and count it in reads
      for (i in 1:length(primersF)) {
        sample <- c(sample, length(grep(paste("^",primersF[i],sep=""), reads, value = TRUE)))
      }
      
      # other
      sample <- c(sample,as.integer(summary(fq)["Length"]) - sum(as.integer(sample[2:length(sample)])))
      
      result="unknown"
      sample <- c(sample,result)
      
      # the first sample defines the data frame
      if (num_sample == 1) {
        forward_table <- data.frame(t(sample))
        colnames(forward_table) <- c("Sample",primersF_names,"Other", "Result")
      }
      # then rows are added
      else{
        forward_table[num_sample,] <- sample
      }
    }
  }
  if (mode=="reverse" | mode=="both") {
    
    for (num_sample in 1:length(reverse_reads)){
      # fastq
      fq <- readFastq(reverse_reads[num_sample])
      
      # sequences
      reads <- sread(fq)
      
      # sample
      sample <- c(substring(basename(reverse_reads[num_sample]), 1, regexpr(reverse_pattern, basename(reverse_reads[num_sample]))-1))
      
      # look at each primer and count it in reads
      for (i in 1:length(primersR)) {
        sample <- c(sample, length(grep(paste("^",primersR[i],sep=""), reads, value = TRUE)))
      }
      
      # other
      sample <- c(sample,as.integer(summary(fq)["Length"]) - sum(as.integer(sample[2:length(sample)])))
      
      result="unknown"
      sample <- c(sample,result)
      
      # the first sample defines the data frame
      if (num_sample == 1) {
        reverse_table <- data.frame(t(sample))
        colnames(reverse_table) <- c("Sample",primersR_names,"Other", "Result")
      }
      # then rows are added
      else{
        reverse_table[num_sample,] <- sample
      }
    }
  }
  
  if (mode=="forward") my_table <- forward_table
  else if (mode=="reverse") my_table <- reverse_table
  else my_table <- list(forward_table,reverse_table)
  return(my_table)
   
}

tables <- cross_contamination(reads_path,primers_path)
if (!(dir.exists((paste(output_path, "primers_check", sep = "/"))))) {
  dir.create(paste(output_path, "primers_check", sep="/"))
}

write.csv(tables[[1]],paste0(output_path,"/primers_check/forward_status.csv"))
write.csv(tables[[2]], paste0(output_path,"/primers_check/reverse_status.csv"))
