#!/bin/bash

source /opt/anaconda3/etc/profile.d/conda.sh

conda activate /home/povp/conda_envs/qiime2-amplicon-2024.5
export TMPDIR=/home/povp/data

path=/home/povp/16S/
path_vysledky=/home/povp/16S_deblur/
path_base=/home/povp/
path_bbmap=/home/povp/bbmap/


cd $path

# removing adapters and phiX
for lib_name in knihovna2 knihovna3 knihovna4 knihovna5 knihovna6 knihovna7 knihovna8 knihovna9
do
        cd ${path}/${lib_name}/
        mkdir ${path_vysledky}/${lib_name}
        mkdir ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk

        for read1 in *R1_001.fastq*
        do
        read2=$(echo $read1| sed 's/R1_/R2_/')
        ${path_bbmap}bbduk.sh \
        ref=${path_bbmap}resources/adapters.fa \
        in1=${read1} in2=${read2} out=stdout.fq \
        k=23 hdist=1 tbo cf=TRUE ftm=5 \
        2> ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/remove_adaptors.log | \
        ${path_bbmap}bbduk.sh \
        ref=${path_bbmap}resources/phix174_ill.ref.fa.gz \
        in=stdin.fq int=t \
        out1=${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/noadaptors_nophix_${read1} \
        out2=${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/noadaptors_nophix_${read2} \
        k=31 hdist=1 2> ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/remove_phiX.log
        done
done
