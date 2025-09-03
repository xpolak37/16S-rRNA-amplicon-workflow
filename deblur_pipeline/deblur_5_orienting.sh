#!/bin/bash

source /opt/anaconda3/etc/profile.d/conda.sh

conda activate /home/povp/conda_envs/qiime2-amplicon-2024.5
export TMPDIR=/home/povp/data

path=/home/povp/16S/
path_vysledky=/home/povp/16S_deblur/
path_base=/home/povp/
path_bbmap=/home/povp/bbmap/


# ORIENTING
conda activate /home/povp/conda_envs/vsearch

for lib_name in knihovna2 knihovna3 knihovna4 knihovna5 knihovna6 knihovna7 knihovna8 knihovna9
do
       cd ${path_vysledky}/${lib_name}/merge_pairs_bbmerge
       mkdir ${path_vysledky}/${lib_name}/oriented
       for read in *mergedpairs.fastq.gz
       do
               sample=$(echo $read| sed 's/trimmed_//')
               sample=$(echo $sample| sed 's/_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz//')
               vsearch --orient $read \
                               --db /home/povp/silva_db/silva_nrr99_v138.1_orienting.udb \
                               --fastqout ${path_vysledky}/${lib_name}/oriented/${sample}_oriented.fq \
                               --tabbedout ${path_vysledky}/${lib_name}/oriented/orient.tx
       done
done
