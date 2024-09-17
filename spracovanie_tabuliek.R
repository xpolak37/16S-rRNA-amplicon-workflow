###############################################################
# change this!!!
path_output="/storage/brno2/home/xpolak37/ikem/spracovane_data/knihovna9/dada_results/"
#############################################################

# ARGS
args = commandArgs(trailingOnly=TRUE)

if (length(args)>0) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

asv_table_path <- paste0(path_output,"final_seqtab.csv")
taxa_table_path <- paste0(path_output,"taxa.csv")

library(biomformat)
library(Biostrings)
library(openxlsx)
library(data.table)
library(dplyr)
library(magrittr)
library(tibble)

# funkcie
build_taxa_table <- function(taxa_table, fasta_file){
  domain <- c()
  phylum <- c()
  class <- c()
  order <- c()
  family <- c()
  genus <- c()
  species <- c()
  for (i in 1:nrow(taxa_table)){
    taxon_string <-taxa_table$Taxon[i]
    taxa <- unlist(strsplit(taxon_string,";"))
    
    domain <- c(domain,ifelse(TRUE %in% grepl("d__",taxa), gsub("d__","",taxa[grep("d__", taxa)]),"NA"))
    phylum <- c(phylum,ifelse(TRUE %in% grepl("p__",taxa), gsub("p__","",taxa[grep("p__", taxa)]),"NA"))
    class <- c(class,ifelse(TRUE %in% grepl("c__",taxa), gsub("c__","",taxa[grep("c__", taxa)]),"NA"))
    order <- c(order,ifelse(TRUE %in% grepl("o__",taxa), gsub("o__","",taxa[grep("o__", taxa)]),"NA"))
    family <- c(family,ifelse(TRUE %in% grepl("f__",taxa), gsub("f__","",taxa[grep("f__", taxa)]),"NA"))
    genus <- c(genus,ifelse(TRUE %in% grepl("g__",taxa), gsub("g__","",taxa[grep("g__", taxa)]),"NA"))
    species <- c(species,ifelse(TRUE %in% grepl("s__",taxa), gsub("s__","",taxa[grep("s__", taxa)]),"NA"))
  }
  
  taxa_table <- data.frame(SeqID=as.character(fasta_file),
                           domain=domain,
                           phylum=phylum,
                           class=class,
                           order=order,
                           family=family,
                           genus=genus,
                           species=species) 
  rownames(taxa_table) <- NULL
  return(taxa_table)
}

merge_taxa <- function(asv_table,taxa_table,taxonomic_level,names=TRUE){
  taxa_ranks <- colnames(taxa_table)[-which(colnames(taxa_table)=="SeqID")]
  where_level <- which(tolower(taxa_ranks)==tolower(taxonomic_level))
  taxa_asv_table <- merge(taxa_table,asv_table, by="SeqID", all=TRUE) 
  taxa_asv_table <- suppressWarnings(taxa_asv_table %>% 
                                       group_by_at(taxa_ranks[1:where_level]) %>%
                                       summarise(across(names(taxa_asv_table)[9:ncol(taxa_asv_table)], sum)))
  
  seq_ids <- apply(taxa_asv_table[,1:where_level],1, function(x){
    a <- paste0(substring(tolower(colnames(taxa_asv_table[,1:where_level])),1,1),"__",x, collapse = ";")
    return(a)
  })
  
  asv_table_sub <- as.data.frame(taxa_asv_table[,(where_level+1):ncol(taxa_asv_table)])
  taxa_table_sub <- as.data.frame(taxa_asv_table[,1:where_level])
  
  asv_table_sub$SeqID <- seq_ids
  taxa_table_sub$SeqID <- seq_ids
  
  asv_table_sub <- asv_table_sub[,c(ncol(asv_table_sub),1:ncol(asv_table_sub)-1)]
  taxa_table_sub <- taxa_table_sub[,c(ncol(taxa_table_sub),1:ncol(taxa_table_sub)-1)]
  return(list(asv_table_sub,taxa_table_sub))
}

spracovanie_tabuliek <- function(asv_table, taxa_table){
  res_species <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Species")
  res_genus <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Genus")
  res_family <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Family")
  res_order <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Order")
  res_class <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Class")
  res_phylum <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Phylum")
  res_domain <- merge_taxa(asv_table,taxa_table,taxonomic_level = "Domain")
  
  # species
  asv_species <- res_species[[1]]
  taxa_species <- res_species[[2]]
  
  # genus
  asv_genus <- res_genus[[1]]
  taxa_genus <- res_genus[[2]]
  
  # family
  asv_family <- res_family[[1]]
  taxa_family <- res_family[[2]]
  
  # order
  asv_order <- res_order[[1]]
  taxa_order <- res_order[[2]]
  
  # class
  asv_class <- res_class[[1]]
  taxa_class <- res_class[[2]]
  
  # phylum
  asv_phylum <- res_phylum[[1]]
  taxa_phylum <- res_phylum[[2]]
  
  # kingdom
  asv_domain <- res_domain[[1]]
  taxa_domain <- res_domain[[2]]
  
  wb <- createWorkbook(); 
  
  addWorksheet(wb, sheetName ="Domain"); 
  writeData(wb, sheet = "Domain", asv_domain)
  setColWidths(wb, c("Domain"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Phylum"); 
  writeData(wb, sheet = "Phylum", asv_phylum)
  setColWidths(wb, c("Phylum"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Class"); 
  writeData(wb, sheet = "Class", asv_class)
  setColWidths(wb, c("Class"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Order"); 
  writeData(wb, sheet = "Order", asv_order)
  setColWidths(wb, c("Order"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Family"); 
  writeData(wb, sheet = "Family", asv_family)
  setColWidths(wb, c("Family"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Genus"); 
  writeData(wb, sheet = "Genus", asv_genus)
  setColWidths(wb, c("Genus"), cols = 1, widths = "auto")
  
  addWorksheet(wb, sheetName ="Species"); 
  writeData(wb, sheet = "Species", asv_species)
  setColWidths(wb, c("Species"), cols = 1, widths = "auto")
  
  return(wb)
}

create_asv_taxa_table <- function(asv_table, taxa_table){
  # creating asv+taxa

  taxa_ranks <- colnames(taxa_table)
  where_level <- which(tolower(taxa_ranks)=="species")
  taxa_asv_table <- merge(taxa_table,asv_table, by="SeqID", all=TRUE) 
  
  seq_ids <- apply(taxa_asv_table[,2:where_level],1, function(x){
    a <- paste0(substring(tolower(colnames(taxa_asv_table[,2:where_level])),1,1),"__",x, collapse = ";")
    return(a)
  })
  
  taxa_asv_table$SeqID <- seq_ids
  taxa_asv_table <- taxa_asv_table[,-which(colnames(taxa_asv_table) %in% taxa_ranks[2:8])]
  return(taxa_asv_table)
}

data_check <- function(asv_table,taxa_table){
  # NaN to 0 and unassigned
  asv_table[is.na(asv_table)] <- 0
  taxa_table[is.na(taxa_table)] <- "unassigned"
  taxa_table[taxa_table ==""] <- "unassigned"
  
  rownames(asv_table) <- NULL
  rownames(taxa_table) <- NULL
  
  # check for redundancy
  to_retain <- rowSums(asv_table[,-1]) != 0
  if(FALSE %in% to_retain){
    to_discard <- asv_table$SeqID[!to_retain]
    message((paste("Removing", sum(!to_retain), "ASV(s)")))
    asv_table <- asv_table[-which(asv_table$SeqID %in% to_discard),]
    #taxa_table <- taxa_table[-which(taxa_table$SeqID %in% to_discard),]
  }
  taxa_table %<>% column_to_rownames("SeqID") 
  taxa_table <- taxa_table[asv_table$SeqID,] %>% rownames_to_column("SeqID")
  
  # removing rownames
  row.names(asv_table) <- NULL
  row.names(taxa_table) <- NULL
  
  return(list(asv_table,taxa_table))
}

## Nacitanie dat na kontrolu ----
asv_tab <- as.data.frame(fread(asv_table_path,header=TRUE))
taxa_table <- as.data.frame(fread(taxa_table_path, header=TRUE))

colnames(asv_tab)[1] <- "SeqID"
colnames(taxa_table)[1] <- "SeqID"

data_checked  <- data_check(asv_table = asv_tab,taxa_table = taxa_table)
asv_tab <- data_checked[[1]]
taxa_table <- data_checked[[2]]

## Vztvorenie tabuliek ----
asv_taxa_table <- create_asv_taxa_table(asv_tab, taxa_table)
wb <-   spracovanie_tabuliek(asv_tab, taxa_table)


write.csv(asv_tab,paste0(path_output,"asv_table.csv"),row.names=FALSE)
write.csv(taxa_table,paste0(path_output,"taxa_table.csv"),row.names=FALSE)
write.csv(asv_taxa_table,paste0(path_output,"asv_taxa_table.csv"), row.names = FALSE)
saveWorkbook(wb,paste0(path_output,"all_levels.xlsx"), overwrite = TRUE)




