#!/bin/bash

source /opt/anaconda3/etc/profile.d/conda.sh

conda activate /home/povp/conda_envs/qiime2-amplicon-2024.5
export TMPDIR=/home/povp/data

path=/home/povp/16S/
path_vysledky=/home/povp/16S_deblur/
path_base=/home/povp/
path_bbmap=/home/povp/bbmap/


# QIIME2

conda activate /home/povp/conda_envs/qiime2-amplicon-2024.5

for lib_name in knihovna2 knihovna3 knihovna4 knihovna5 knihovna6 knihovna7 knihovna8
do
	cd ${path_vysledky}/${lib_name}/oriented
	mkdir ../qiime2
	mkdir ../qiime2/output
	echo -e "sample-id\tabsolute-filepath\tdirection" > ../qiime2/manifest_${lib_name}.tsv
	
	for file in *_*_oriented.fq; do 
		mv "$file" "${file//_/-}"
	done
	
	
	for read in *-oriented.fq; do
		sample=$(echo $read| sed 's/-oriented.fq//')
		echo -e "${sample}\t${path_vysledky}${lib_name}/oriented/${read}\tforward" >>  ../qiime2/manifest_${lib_name}.tsv
	done

	cd ../qiime2/

	qiime tools import \
    	--input-path manifest_${lib_name}.tsv \
    	--type 'SampleData[SequencesWithQuality]' \
    	--input-format SingleEndFastqManifestPhred33V2 \
    	--output-path ${lib_name}.qza

    	qiime demux summarize --i-data ${lib_name}.qza --o-visualization output/${lib_name}_imported_seqs.qzv --p-n 10000

    	qiime deblur denoise-16S \
	--i-demultiplexed-seqs ${lib_name}.qza \
	--p-trim-length 400 \
	--p-sample-stats \
	--o-representative-sequences output/${lib_name}_ASV_seqs.qza \
	--o-table output/${lib_name}_ASV_abundance.qza \
	--o-stats output/${lib_name}_stats_deblur.qza
	
	qiime feature-classifier classify-sklearn \
	--i-reads ${path_vysledky}/${lib_name}/qiime2/output/${lib_name}_ASV_seqs.qza \
	--i-classifier /home/povp/Rescript_classifier/classifier_138_1_new/silva-138.1-ssu-nr99-319f-806r-classifier.qza \
	--p-n-jobs 40 \
	--o-classification ${path_vysledky}/${lib_name}/qiime2/output/${lib_name}_ASV_taxonomy.qza
	done


mkdir /home/povp/16S_deblur/vysledky

qiime feature-table merge \
    --i-tables /home/povp/16S_deblur/knihovna2/qiime2/output/knihovna2_ASV_abundance.qza \
    --i-tables /home/povp/16S_deblur/knihovna3/qiime2/output/knihovna3_ASV_abundance.qza \
    --i-tables /home/povp/16S_deblur/knihovna4/qiime2/output/knihovna4_ASV_abundance.qza \
    --i-tables /home/povp/16S_deblur/knihovna5/qiime2/output/knihovna5_ASV_abundance.qza \
    --i-tables /home/povp/16S_deblur/knihovna6/qiime2/output/knihovna6_ASV_abundance.qza \
    --i-tables /home/povp/16S_deblur/knihovna7/qiime2/output/knihovna7_ASV_abundance.qza \
    --i-tables /home/povp/16S_deblur/knihovna8/qiime2/output/knihovna8_ASV_abundance.qza \
    --i-tables /home/povp/16S_deblur/knihovna9/qiime2/output/knihovna9_ASV_abundance.qza \
    --o-merged-table /home/povp/16S_deblur/vysledky/ASV_abundance_merged.qza

qiime feature-table merge-seqs \
    --i-data /home/povp/16S_deblur/knihovna2/qiime2/output/knihovna2_ASV_seqs.qza \
    --i-data /home/povp/16S_deblur/knihovna3/qiime2/output/knihovna3_ASV_seqs.qza \
    --i-data /home/povp/16S_deblur/knihovna4/qiime2/output/knihovna4_ASV_seqs.qza \
    --i-data /home/povp/16S_deblur/knihovna5/qiime2/output/knihovna5_ASV_seqs.qza \
    --i-data /home/povp/16S_deblur/knihovna6/qiime2/output/knihovna6_ASV_seqs.qza \
    --i-data /home/povp/16S_deblur/knihovna7/qiime2/output/knihovna7_ASV_seqs.qza \
    --i-data /home/povp/16S_deblur/knihovna8/qiime2/output/knihovna8_ASV_seqs.qza \
    --i-data /home/povp/16S_deblur/knihovna9/qiime2/output/knihovna9_ASV_seqs.qza \
    --o-merged-data /home/povp/16S_deblur/vysledky/ASV_seqs_merged.qza

qiime feature-table merge-taxa \
    --i-data /home/povp/16S_deblur/knihovna2/qiime2/output/knihovna2_ASV_taxonomy.qza \
    --i-data /home/povp/16S_deblur/knihovna3/qiime2/output/knihovna3_ASV_taxonomy.qza \
    --i-data /home/povp/16S_deblur/knihovna4/qiime2/output/knihovna4_ASV_taxonomy.qza \
    --i-data /home/povp/16S_deblur/knihovna5/qiime2/output/knihovna5_ASV_taxonomy.qza \
    --i-data /home/povp/16S_deblur/knihovna6/qiime2/output/knihovna6_ASV_taxonomy.qza \
    --i-data /home/povp/16S_deblur/knihovna7/qiime2/output/knihovna7_ASV_taxonomy.qza \
    --i-data /home/povp/16S_deblur/knihovna8/qiime2/output/knihovna8_ASV_taxonomy.qza \
    --i-data /home/povp/16S_deblur/knihovna9/qiime2/output/knihovna9_ASV_taxonomy.qza \
    --o-merged-data /home/povp/16S_deblur/vysledky/ASV_taxonomy_merged.qza

cd /home/povp/16S_deblur/vysledky

qiime taxa filter-table \
	--i-table ASV_abundance_merged.qza \
	--i-taxonomy ASV_taxonomy_merged.qza \
	--p-include 'o__' \
	--p-mode 'contains' \
	--p-exclude 'mitochondria,chloroplast,o__;' \
	--o-filtered-table ASV_abundance_filtered.qza

mkdir /home/povp/16S_deblur/vysledky/exported

qiime tools export \
--input-path ASV_abundance_filtered.qza \
--output-path exported
qiime tools export \
--input-path ASV_seqs_merged.qza \
--output-path exported
qiime tools export \
--input-path ASV_taxonomy_merged.qza \
--output-path exported
