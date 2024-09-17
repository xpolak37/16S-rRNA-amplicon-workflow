# 16S rRNA amplicon workflow

This repository includes whole bioinformatics workflow for 16S rRNA amplicon reads processing,
all parameters are set to work properly on specific data from our lab. Before any usage, please
check if it fits your data as well. 

Whole workflow can be execute main.sh, where you have to change paths to conda environments and to your data
```bash
bash main.sh
```

**1. Cross contamination**

cross_contamination.R creates a report in form of data.frame with individual primers counts in reads across all samples

**2. Cutadapt**

cutadapt.sh script checks the quality of the raw reads and creates the MULTIQC report, run cutadapt for adapter and primer trimming. Lastly,
it creates MULTIQC report for trimmed reads as well. 

**3. Merging + orienting**

merging_orienting.sh script use VSEARCH to merge the paired-end reads and to orient them to have the same direction for further processing.

**4. DADA2**

dada_merged_reads.R script includes DADA2 pipeline for merged reads (single-end) as well as the taxonomic classification using IDTAXA and addspecies()

**5. Results processing**

results_processing.R script processes the results of DADA2 pipeline into multiple dataframes and excel sheets:
- asv_table.csv: table with counts for each ASV
- taxa_table.csv: table with taxonomy for each ASV
- asv_taxa_table.csv: table with counts for each taxonomy
- all_levels.xlsx: excel sheets with counts aggregated to each taxonomic level
