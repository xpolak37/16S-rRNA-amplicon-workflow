#!/bin/bash
eval "$(conda shell.bash hook)"

# variables for conda environment
conda_env_dir_preprocessing="/home/povp/conda_envs/WGS_preprocessing/R_v4_env"
conda_env_dir_R="/home/povp/conda_envs/R_v4_env"

# paths for data storage
path_input="/home/povp/seq_data/16S/knihovna11/decontaminated/human_phix"
path_output="/home/povp/16S_deblur/knihovna11/dada2/"
path_project_dir="/home/povp/16S_deblur/knihovna11"

# activate environment
conda activate ${conda_env_dir_preprocessing}

# info for tools and versions txt
echo "main_dada2.sh:" >> ${path_project_dir}/run_info/tools.txt

Rscript dada2_processing.R ${path_input} ${path_output} ${path_project_dir}

