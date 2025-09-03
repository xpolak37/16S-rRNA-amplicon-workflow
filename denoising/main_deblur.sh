#!/bin/bash
eval "$(conda shell.bash hook)"

# variables for conda environment
conda_env_dir_preprocessing="/home/povp/conda_envs/WGS_preprocessing/"
conda_env_dir_R="/home/povp/conda_envs/R_v4_env"
conda_env_dir_qiime="/home/povp/conda_envs/qiime2-amplicon-2024.5"

# paths for data storage
export path_input="/home/povp/seq_data/16S/knihovna11/decontaminated/human_phix"
export path_output="/home/povp/16S_deblur/knihovna11/deblur/"
export path_project_dir="/home/povp/16S_deblur/knihovna11"/
export path_bbmap=/home/povp/progs/bbmap/
export TMPDIR=/home/povp/tmp

# activate environment
conda activate ${conda_env_dir_preprocessing}

# info for tools and versions txt
echo "main_deblur.sh:" >> ${path_project_dir}/run_info/tools.txt

# making dirs
mkdir ${path_output}
mkdir ${path_output}/merged
mkdir ${path_output}/oriented
mkdir ${path_output}/qiime2
mkdir ${path_output}/qiime2/output
path_qiime_input=${path_output}/merged

# MERGING READS
run_bbmap(){
        read1=$1
        read2=$(echo $read1| sed 's/R1/R2/')
        sample=$(echo "$read1" | sed 's/_L[0-9][0-9][0-9].*//')

        ${path_bbmap}bbmerge.sh \
        in1=${read1} \
        in2=${read2} \
        out=stdout.fq qtrim=r trimq=15 maxlength=500 mininsert=350 | \
        ${path_bbmap}reformat.sh -Xmx20g -int=f -maxns=0 \
        in=stdin.fq out=${path_output}/merged/${sample}_mergedpairs.fastq.gz
        
}

find ${path_input} -type f -name "*R1_trimmed_cleaned.fastq.gz" | parallel -j 10 run_bbmap


# ORIENTING READS ?????????

# DEBLUR QIIME2
conda activate ${conda_env_dir_qiime}

# samplesheet
echo -e "sample-id\tabsolute-filepath\tdirection" > ${path_output}/qiime2/manifest.tsv
	
for file in ${path_qiime_input}/*_*_oriented.fq; do 
        mv "$file" "${file//_/-}"
done
	
	
for read in ${path_qiime_input}/*-oriented.fq; do
	sample=$(echo $read| sed 's/-oriented.fq//')
	echo -e "${sample}\t${path_output}/oriented/${read}\tforward" >>  ${path_output}/qiime2/manifest.tsv
done

cd ${path_output}/qiime2

# Import
qiime tools import \
        --input-path manifest.tsv \
    	--type 'SampleData[SequencesWithQuality]' \
    	--input-format SingleEndFastqManifestPhred33V2 \
    	--output-path imported_data.qza

# Deblur
qiime deblur denoise-16S \
        --i-demultiplexed-seqs imported_data.qza \
	--p-trim-length 400 \
	--p-sample-stats \
	--o-representative-sequences output/${lib_name}_ASV_seqs.qza \
	--o-table output/ASV_abundance.qza \
	--o-stats output/stats_deblur.qza
	

#qiime taxa filter-table \
#	--i-table ASV_abundance.qza \
#	--i-taxonomy ASV_taxonomy_merged.qza \
#	--p-include 'o__' \
#	--p-mode 'contains' \
#	--p-exclude 'mitochondria,chloroplast,o__;' \
#	--o-filtered-table ASV_abundance_filtered.qza

qiime tools export \
--input-path output/ASV_abundance.qza \
--output-path exported
qiime tools export \
--input-path output/ASV_seqs_merged.qza \
--output-path exported
qiime tools export \
--input-path ASV_taxonomy_merged.qza \
--output-path exported
