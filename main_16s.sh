#!/bin/bash
set -euo pipefail

eval "$(conda shell.bash hook)"

perform_quality_control=TRUE
perform_preprocessing=TRUE
perform_denoising=TRUE
perform_taxassignment=TRUE
perform_dada2=FALSE
perform_deblur=FALSE
keep_intermediate=FALSE

# DIRECTORIES
input_dir="/home/povp/seq_data/16S/knihovna11/data/"
project_dir="/home/povp/16S_deblur/knihovna11"

# CONDA ENVIRONMENTS
quality_env="/home/povp/conda_envs/quality"
blast_env="/home/povp/conda_envs/blast_env"
R_env="/home/povp/conda_envs/R_v4_env"
preprocessing_env="/home/povp/conda_envs/WGS_preprocessing"
nextflow_env="/home/povp/conda_envs/WGS_preprocessing"
metaphlan_env="/home/povp/conda_envs/WGS"

# set all variables
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir) input_dir=$2; shift 2 ;;
        --project_dir) project_dir=$2; shift 2 ;;
        --keep_intermediate) keep_intermediate=TRUE; shift ;;
        --skip_quality) perform_quality_control=FALSE; shift ;;
        --skip_preprocessing) perform_preprocessing=FALSE; shift ;;
        --skip_taxassignment) perform_taxprofiling=FALSE; shift ;;
        --dada2) perform_dada2=TRUE; shift ;;
        --deblur) perform_deblur=TRUE; shift ;;
        --quality_env) quality_env=$2; shift 2;;
        --blast_env) blast_env=$2; shift 2;;
        --R_env) R_env=$2; shift 2;;
        --preprocessing_env) preprocessing_env=$2; shift 2;;
        --qiime_env) metaphlan_env=$2; shift 2;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

export path_scripts="$(pwd)"

# define the txt file for version tracking
mkdir -p ${project_dir}/run_info/
echo "# Tools and their versions throughout the run" > ${project_dir}/run_info/tools.txt

# RUN the pipeline

## Quality control
if [[ "${perform_quality_control}" == "TRUE" ]]; then
    echo "Performing quality control..."
    cd quality_control
    bash main_qc.sh
    echo "Performing cross contamination check"
    Rscript cross_contamination.R ${input_dir} ${path_scripts}/primers.fasta ${project_dir}
    echo "Done"
    cd ${path_scripts}
fi

# Preprocessing
if [[ "${perform_preprocessing}" == "TRUE" ]]; then
    cd preprocessing
    echo "Performing trimming with cutadapt"
    bash cutadapt.sh
    echo "Performing host decontamination with hostile"
    bash host_decontamination.sh
    echo "Done"
    cd ${path_scripts}
fi

# Denoising
if [[ "${perform_denoising}" == "TRUE" ]]; then
    if [[ "$perform_dada2" == "FALSE" && "$perform_deblur" == "FALSE" ]]; then
        perform_dada2=TRUE
    fi

    cd denoising

    if [[ "${perform_dada2}" == "TRUE" ]]; then
        echo "Running Dada2 for denoising ..."
        bash main_dada2.sh
        echo Done
    fi
    if [[ "${perform_deblur}" == "TRUE" ]]; then
        echo "Running Deblur for denoising ..."
        bash main_deblur.sh
        echo Done
    fi
    cd ${path_scripts}
fi

# Taxonomic assignment
if [[ "${perform_taxassignment}" == "TRUE" ]]; then
    if [[ "$perform_assigntaxa" == "FALSE" && "$perform_IDTAXA" == "FALSE" ]]; then
        perform_IDTAXA=TRUE
    fi

    cd taxonomic_assignment
    echo "Running taxonomic assignment"
    bash main_tax.sh
    echo Done
    if [[ "${perform_assigntaxa}" == "TRUE" ]]; then
        echo "Running Dada2 for denoising ..."
        bash main_dada2.sh
        echo Done
    fi
    if [[ "${perform_IDTAXA}" == "TRUE" ]]; then
        echo "Running Deblur for denoising ..."
        bash main_deblur.sh
        echo Done
    fi
    cd ${path_scripts}
fi


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
