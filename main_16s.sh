#!/bin/bash

eval "$(conda shell.bash hook)"

####### USER INPUTS #########
# DIRECTORIES
export input_dir="/home/povp/seq_data/16S/knihovna8/renaming/"
export project_dir="/home/povp/16S_deblur/knihovna8_t1d/"
export path_bbmap="/home/povp/progs/bbmap/"
export TMPDIR="/home/povp/tmp"
export orienting_db="/home/povp/taxonomic_classifiers/vsearch_classifier/silva_nrr99_v138.1_orienting.udb"
export qiime_classifier="/home/povp/taxonomic_classifiers/Rescript_classifier/classifier_341f_805r/silva-138.1-ssu-nr99-341f-805r-classifier.qza"
export decipher_classifier="/home/povp/taxonomic_classifiers/decipher_classifier/idtaxa_trainingSet_V3V4_silva_138_2.RData"

# CONDA ENVIRONMENTS
export amplicon_16s_env="/home/povp/conda_envs/16S_amplicon/"
export blast_env="/home/povp/conda_envs/blast_env"
export qiime_env="/home/povp/conda_envs/qiime2-amplicon-2024.5"

# PARAMETERS
## BBMAP merging and filtering
export maxns=1
export trimq=15
export maxlength=500
export mininsert=350

## DEBLUR
export p_trim_length=400
export min_reads=1
export min_size=1

## IDTAXA
export threshold=60
export strand="both"

##############################

conda activate ${amplicon_16s_env}

perform_quality_control=TRUE
perform_preprocessing=TRUE
perform_denoising=TRUE
perform_taxassignment=TRUE
perform_dada2=FALSE
perform_deblur=FALSE
keep_intermediate=FALSE
perform_assigntaxa=FALSE
perform_IDTAXA=FALSE
perform_qiime_tax=FALSE

# set all variables
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir) input_dir=$2; shift 2 ;;
        --project_dir) project_dir=$2; shift 2 ;;
        --keep_intermediate) keep_intermediate=TRUE; shift ;;
        --skip_quality) perform_quality_control=FALSE; shift ;;
        --skip_preprocessing) perform_preprocessing=FALSE; shift ;;
        --skip_denoising) perform_denoising=FALSE; shift ;;
        --skip_taxassignment) perform_taxassignment=FALSE; shift ;;
        --dada2) perform_dada2=TRUE; shift ;;
        --deblur) perform_deblur=TRUE; shift ;;
        --qiime_tax) perform_qiime_tax=TRUE; shift ;;
        --idtaxa) perform_IDTAXA=TRUE; shift;;
        --assigntaxa) perform_assigntaxa=TRUE; shift;;
        --blast_env) blast_env=$2; shift 2;;
        --main_env) amplicon_16s_env=$2; shift 2;;
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
    bash main_qc.sh ${amplicon_16s_env} ${blast_env} ${input_dir} ${project_dir}
    conda activate ${amplicon_16s_env}
    echo "Performing cross contamination check"
    Rscript cross_contamination.R ${input_dir} ${path_scripts}/primers.fasta ${project_dir}
    echo "Done"
    cd ${path_scripts}
fi

# Preprocessing
if [[ "${perform_preprocessing}" == "TRUE" ]]; then
    cd preprocessing
    echo "Performing trimming with cutadapt"
    bash cutadapt.sh ${amplicon_16s_env} ${input_dir} ${project_dir}
    echo "Performing host decontamination with hostile"
    bash host_decontamination.sh ${amplicon_16s_env} ${project_dir}/preprocessed/ ${project_dir}
    echo "Done"
    cd ${path_scripts}
fi

# Denoising
if [[ "${perform_denoising}" == "TRUE" ]]; then
    if [[ "$perform_dada2" == "FALSE" && "$perform_deblur" == "FALSE" ]]; then
        perform_dada2=TRUE
    fi

    cd denoising

    mkdir ${project_dir}/denoised
    if [[ "${perform_dada2}" == "TRUE" ]]; then
        echo "Running Dada2 for denoising ..."
        bash main_dada2.sh ${project_dir}/decontaminated/ ${project_dir}
        echo Done
        tax_input=${project_dir}/denoised/dada2
    fi
    if [[ "${perform_deblur}" == "TRUE" ]]; then
        echo "Running Deblur for denoising ..."
        bash main_deblur.sh ${project_dir}/decontaminated/ ${project_dir}
        echo Done
        tax_input=${project_dir}/denoised/deblur
    fi
    cd ${path_scripts}
fi

# Taxonomic assignment
if [[ "${perform_taxassignment}" == "TRUE" ]]; then

    # make output directory
    mkdir ${project_dir}/taxonomy

    # when no classifier preference was set, perform idtaxa()
    if [[ "${perform_assigntaxa}" == "FALSE" && "${perform_IDTAXA}" == "FALSE" && "${perform_qiime_tax}" == "FALSE" ]]; then
        perform_IDTAXA=TRUE
    fi

    # if there was no denoising step, assume the deblur/dada2 was previously run
    if [[ "${perform_denoising}" == "FALSE" ]]; then

        # use deblur output as default
        tax_input=${project_dir}/denoised/deblur

        # change according to preference
        if [[ "${perform_dada2}" == "TRUE" ]]; then
            tax_input=${project_dir}/denoised/dada2
        fi
        if [[ "${perform_deblur}" == "TRUE" ]]; then
            tax_input=${project_dir}/denoised/deblur
        fi
    fi

    cd taxonomic_assignment
    if [[ "${perform_assigntaxa}" == "TRUE" ]]; then
        echo "Running assignTaxa() for taxonomic assignment ..."

        # make directory
        mkdir ${project_dir}/taxonomy/assigntaxa

        # run script
        Rscript assign_taxonomy.R
        echo Done
    fi
    if [[ "${perform_IDTAXA}" == "TRUE" ]]; then
        echo "Running idtaxa() for taxonomic assignment ..."

        # make dir
        mkdir ${project_dir}/taxonomy/idtaxa

        # run script
        Rscript decipher_idtaxa.R ${tax_input} ${project_dir}/taxonomy/idtaxa ${project_dir} ${decipher_classifier}
        echo Done

        # cp the result in run_info
        cp ${project_dir}/taxonomy/idtaxa/taxa_table.csv ${project_dir}/run_info/taxa_table.csv
        cp ${project_dir}/taxonomy/idtaxa/taxa_table_conf.csv ${project_dir}/run_info/taxa_table_conf.csv

    fi
    cd ${path_scripts}
fi


#echo "Starting step1 - cross_contamination.R"
# 1st step: statistics for primer pairs
#Rscript cross_contamination.R "reads_path='$reads_path'" "primers_path='$primers_path'" "output_path='$output_path'"

#echo "Done"

#echo "Starting step2 - cutadapt.sh"
# 2nd step: cutadapt + quality score
#bash cutadapt.sh --path $reads_path --path_output $output_path --path_base $path_base

#echo "Done"

#echo "Starting step3- merging_orienting.sh"
# 3th step: merging and orienting
#bash merging_orienting.sh --path $reads_path --path_output $output_path --path_base $path_base

#echo "Done"

#echo "Starting step4 - dada_merged_reads.R"
# 4th step: DADA
#Rscript dada_merged_reads.R "path='$output_path/merged_oriented/'" "path_to_results='$output_path/dada_results/'" "path_idtaxa='$path_idtaxa'" "path_addspecies='$path_addspecies'"

#echo "Done"

#echo "Starting step5 - results_processing.R"
# 5th step: creating dataframes
#Rscript results_processing.R "path_output='$output_path/dada_results/'"

#echo "Juhu - koniec"
