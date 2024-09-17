#!/bin/bash
conda init --all
source ~/.bashrc

# conda environments - change this!
#################################################################
conda_envs_dir=/home/povp/conda_envs
r_env=/home/povp/conda_evs/R_v4_env
vsearch_env=/home/povp/conda_evs/vsearch
###################################################################

###############################################################
# change this!!!!
path="home/povp/16S/knihovna9/"
path_vysledky="home/povp/spracovane_data/knihovna9/"
path_base="home/povp/"

#############################################################

# arguments

while [[ $# -gt 0 ]]; do
 case $1 in
    --path)
       path_arg="$2"
       shift
       shift
       ;;
    --path_output)
       path_output_arg="$2"
       shift
       shift
       ;;
    --path_base)
       bath_base_arg="$2"
       shift
       shift
       ;;
     *)
    echo "Input not understood."
       exit 1
       ;;
 esac
done

echo "Starting cutadapt.sh with arguments:"
# arguments control
if [ -z ${path_arg+x} ]; then path_arg=$path; fi
if [ -z ${path_output_arg+x} ]; then path_output_arg=$path_output; fi 
if [ -z ${path_base_arg+x} ]; then path_base_arg=$path_base; fi 

echo "Path: ${path_arg}"
echo "Path output: ${path_output_arg}"
echo "Path base: ${path_base_arg}"

path=$path_arg
path_output=$path_output_arg
path_base=$path_base_arg

# merging and orienting reads
cd $conda_envs_dir
conda activate $vsearch_env

cd $path_output/trimmed
mkdir $path_output/merged
mkdir $path_output/merged_oriented

for read1 in *R1_001.fastq*
do
read2=$(echo $read1| sed 's/R1_/R2_/')
sample=$(basename "$read1" | cut -d '_' -f 2,3)

vsearch --fastq_mergepairs ${read1} --reverse ${read2} --fastqout ../merged/${sample}_paired.fastq 

vsearch --orient ../merged/${sample}_paired.fastq --db /home/povp/silva_db/silva_nrr99_v138.1_orienting.udb --fastqout ../merged_oriented/${sample}_oriented.fq 
done
