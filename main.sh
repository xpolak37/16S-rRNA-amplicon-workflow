#!/bin/bash
source activate base
conda init
source ~/.bashrc

# bash script for processing 16S rRNA amplicon sequencing data

# conda environments
#################################################################
r_env=/path/to/conda_environment
###################################################################

# Main paths
################################################################################
scripts_path=/path/to/this/repo/
reads_path=/path/to/reads/directory
primers_path=/path/to/FASTA/primers
output_path=/path/to/output/directory
path_idtaxa=/path/to/SILVA_SSU_r138_2019.RData
path_addspecies=/path/to/silva_species_assignment_v138.1.fa.gz
################################################################################

conda activate r_env 
cd $scripts_path

echo "Starting step1 - cross_contamination.R"
# 1st step: statistics for primer pairs
Rscript cross_contamination.R "reads_path='$reads_path'" "primers_path='$primers_path'" "output_path='$output_path'"

echo "Done"

echo "Starting step2 - cutadapt.sh"
# 2nd step: cutadapt + quality score
bash cutadapt.sh --path $reads_path --path_output $output_path --path_base $path_base

echo "Done"

echo "Starting step3- merging_orienting.sh"
# 3th step: merging and orienting
bash merging_orienting.sh --path $reads_path --path_output $output_path --path_base $path_base

echo "Done"

echo "Starting step4 - dada_merged_reads.R"
# 4th step: DADA
Rscript dada_merged_reads.R "path='$output_path/merged_oriented/'" "path_to_results='$output_path/dada_results/'" "path_idtaxa='$path_idtaxa'" "path_addspecies='$path_addspecies'"

echo "Done"

echo "Starting step5 - results_processing.R"
# 5th step: creating dataframes
Rscript results_processing.R "path_output='$output_path/dada_results/'"

echo "Juhu - koniec"
