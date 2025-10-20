#!/bin/bash
eval "$(conda shell.bash hook)"

# paths for data storage
export path_input=$1
export path_project_dir=$2
export path_output=${path_project_dir}/denoised/

# activate environment
conda activate ${amplicon_16s_env}

# info for tools and versions txt
echo -e "\nmain_deblur.sh:" >> ${path_project_dir}/run_info/tools.txt

# making dirs
mkdir ${path_output}
mkdir ${path_output}/merged
mkdir ${path_output}/oriented
mkdir ${path_output}/deblur/
mkdir ${path_output}/deblur/qiime2
mkdir ${path_output}/deblur/qiime2/output

# Preprocessing
bash qiime_preprocessing.sh ${path_input} ${path_project_dir} ${path_output}

# QIIME + DEBLUR
path_qiime_input=${path_output}/oriented

path_output=${path_output}/deblur

# DEBLUR QIIME2
conda activate ${qiime_env}

cd ${path_qiime_input}

for f in *_*; do
  [ -f "$f" ] && mv "$f" "${f//_/-}"
done

# samplesheet
echo -e "sample-id\tabsolute-filepath\tdirection" > ${path_output}/qiime2/manifest.tsv
	
for read in ${path_qiime_input}/*-oriented.fq; do
	sample=$(echo $read| sed 's/-oriented.fq//')
        sample=$(basename ${sample})
	echo -e "${sample}\t${read}\tforward" >>  ${path_output}/qiime2/manifest.tsv
done

# QIIME


## Dereplicate sequences to get FeatureData[Sequence] which will be used as reference
## db in deblur to avoid positive filtering. 
#cd ${path_scripts}/denoising
#python build_reference_db.py ${path_project_dir}/oriented/ ${path_project_dir}/qiime2/reference_seqs.fasta

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
	--p-trim-length ${p_trim_length} \
        --p-min-reads ${min_reads} \
        --p-min-size ${min_size} \
	--p-sample-stats \
	--o-representative-sequences output/ASV_seqs.qza \
	--o-table output/ASV_abundance.qza \
	--o-stats output/stats_deblur.qza \
        --p-jobs-to-start 20

qiime tools export \
--input-path output/ASV_abundance.qza \
--output-path output/exported
qiime tools export \
--input-path output/ASV_seqs.qza \
--output-path output/exported
qiime tools export \
--input-path output/stats_deblur.qza \
--output-path output/exported

biom convert -i output/exported/feature-table.biom -o output/exported/feature-table.tsv --to-tsv

## track version
qiime --version >> ${path_project_dir}/run_info/tools.txt

# POSTPROCESSING
conda activate ${amplicon_16s_env}
cd ${path_scripts}/denoising
Rscript qiime_postprocessing.R ${path_output}/qiime2/output/exported ${path_output}

## track version
R --version | head -n 1 >> ${path_project_dir}/run_info/tools.txt

# final counts
## Output file
outfile=${path_project_dir}/run_info/read_counts_summary.txt
helper_outfile=${path_project_dir}/run_info/seqkit_output.txt

## Header
echo -e "Sample\tfastqc-total_sequences" > ${helper_outfile}
awk 'BEGIN{OFS="\t"} NR==1{$(NF+1)="Final_deblur"} NR>1{$(NF+1)=""}1' ${outfile} > tmp && mv tmp ${outfile}

## read count
awk '
NR==1 {
    for (i=2; i<=NF; i++) samples[i]=$i
    next
}
{
    for (i=2; i<=NF; i++) depth[i]+=$i
}
END {
    for (i=2; i<=NF; i++) print samples[i] "\t" depth[i]
}' "${path_output}/asv_table.tsv" >> ${helper_outfile}

cd ${path_scripts}/preprocessing
python3 counts_summary.py ${helper_outfile} ${outfile} 5

## track version
python --version >> ${path_project_dir}/run_info/tools.txt

## remove helper file
rm ${helper_outfile}

# cp the result in run_info
cp ${path_output}/asv_table.tsv ${path_project_dir}/run_info/asv_table.tsv
