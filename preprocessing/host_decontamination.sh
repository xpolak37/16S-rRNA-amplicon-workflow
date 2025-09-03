#!/bin/bash
eval "$(conda shell.bash hook)"

keep_intermediate=FALSE

# variables for conda environment
conda_env_dir_preprocessing="/home/povp/conda_envs/WGS_preprocessing"

# paths for data storage
path_input="/home/povp/16S_deblur/knihovna11/trimmed/"
path_output="/home/povp/16S_deblur/knihovna11/decontaminated/"
path_project_dir="/home/povp/16S_deblur/knihovna11"

# path to phix bowtie indexed
path_bowtie_phix="/home/povp/scripts/Metagenomics_workflow/preprocessing/bowtie2_indexes/phiX174"

# activate environment
conda activate ${conda_env_dir_preprocessing}

# TMPDIR
export TMPDIR=/home/povp/tmp

# info for tools and versions txt
echo "host_decontamination.sh:" >> ${path_project_dir}/run_info/tools.txt

# HOST DECONTAMINATION WILL INCLUDE TWO STEPS:
# 1. Human decontamination using HOSTILE
# 2. PHIX decontamination using HOSTILE

mkdir ${path_output}
mkdir ${path_output}/human
mkdir ${path_output}/human_phix

# human decontamination
find ${path_input} -type f -name "*_R1_trimmed.fastq.gz" | sed 's/_R1_trimmed.fastq.gz//' | parallel -j 10 "hostile clean \
--fastq1 {}_R1_trimmed.fastq.gz \
--fastq2 {}_R2_trimmed.fastq.gz \
--index human-t2t-hla-argos985-mycob140 \
--output ${path_output}/human/ \
--threads 5"

# REMOVE FASTP TRIMMING RESULTS
if [[ "${keep_intermediate}" == "FALSE" ]]; then
    rm -r ${path_project_dir}/trimmed
fi

# phix decontamination
find ${path_output}/human/ -type f -name "*_R1_trimmed.clean_1.fastq.gz" | sed 's/_R1_trimmed.clean_1.fastq.gz//' | parallel -j 10 "hostile clean \
--fastq1 {}_R1_trimmed.clean_1.fastq.gz \
--fastq2 {}_R2_trimmed.clean_2.fastq.gz \
--index ${path_bowtie_phix} \
--output ${path_output}/human_phix/ \
--threads 5"

## track version
echo "hostile" $(hostile --version) >> ${path_project_dir}/run_info/tools.txt

# REMOVE HUMAN DECONTAMINATION RESULTS
if [[ "${keep_intermediate}" == "FALSE" ]]; then
    rm -r ${path_project_dir}/decontaminated/human/
fi

# renaming samples to end with '_trimmed_cleaned.fastq.gz'
for f in ${path_output}/human_phix/*_trimmed.clean_*.fastq.gz; do
    newname=$(echo "$f" | sed -E 's/_trimmed\.clean_[12](\.clean_[12])?/_trimmed_cleaned/')
    mv "$f" "$newname"
done

# The final counts statistics
## Output file
outfile=${path_project_dir}/run_info/read_counts_summary.txt
helper_outfile=${path_project_dir}/run_info/seqkit_output.txt

# Header
echo -e "Sample\tfastqc-total_sequences" > ${helper_outfile}

# read count
ls ${path_output}/human_phix/*.fastq.gz | parallel -j 8 'echo -e "$(basename {})\t$(gzip -cd {} | wc -l | awk "{print \$1/4}")"' >> ${helper_outfile}

cd ${path_scripts}/preprocessing
python3 counts_summary.py ${helper_outfile} ${outfile} 3
## track version
python --version >> ${path_project_dir}/run_info/tools.txt

# remove helper file
rm ${helper_outfile}