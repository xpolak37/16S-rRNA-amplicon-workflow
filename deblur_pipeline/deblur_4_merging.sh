#!/bin/bash

source /opt/anaconda3/etc/profile.d/conda.sh

conda activate /home/povp/conda_envs/qiime2-amplicon-2024.5
export TMPDIR=/home/povp/data

path=/home/povp/16S/
path_vysledky=/home/povp/16S_deblur/
path_base=/home/povp/
path_bbmap=/home/povp/bbmap/


# MERGING
for lib_name in knihovna2 knihovna3 knihovna4 knihovna5 knihovna6 knihovna7 knihovna8 knihovna9
do
       cd ${path_vysledky}/${lib_name}/trimmed
       mkdir ${path_vysledky}/${lib_name}/merge_pairs_bbmerge

       for read1 in *R1_001.fastq*
       do

               read2=$(echo $read1| sed 's/R1/R2/')
               sample1=$(echo $read1| sed 's/trimmed_//')
               sample2=$(echo $read2| sed 's/trimmed_//')
               sample=$(echo $read1| sed 's/_L001_R1_001.fastq.gz//')
               ${path_bbmap}bbmerge.sh \
               in1=${path_vysledky}/${lib_name}/trimmed/${read1} \
               in2=${path_vysledky}/${lib_name}/trimmed/${read2} \
       out=stdout.fq qtrim=r trimq=15 maxlength=500 mininsert=350 ihist=../merge_pairs_bbmerge/${sample}_insert_histogram.ihist \
       | ${path_bbmap}reformat.sh -Xmx20g -int=f -maxns=0 \
       in=stdin.fq out=${path_vysledky}/${lib_name}/merge_pairs_bbmerge/${sample}_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz \
       lhist=${path_vysledky}/${lib_name}/merge_pairs_bbmerge/${sample}_histogram_after_N_filter.txt;
       done
done
