#!/bin/bash
eval "$(conda shell.bash hook)"

# paths for data storage
export path_input=$1
export path_project_dir=$2
export path_output=$3

# activate environment
conda activate ${amplicon_16s_env}

# MERGING READS
run_bbmap(){
        read1=$1
        read2=$(echo $read1| sed 's/R1/R2/')
        sample=$(echo "$read1" | sed 's/_L[0-9][0-9][0-9].*//')
        sample=$(basename ${sample})
        ${path_bbmap}bbmerge.sh \
        in1=${read1} \
        in2=${read2} \
        out=stdout.fq qtrim=r trimq=${trimq} maxlength=${maxlength} mininsert=${mininsert} | \
        ${path_bbmap}reformat.sh -Xmx20g -int=f -maxns=${maxns} \
        in=stdin.fq out=${path_output}/merged/${sample}_mergedpairs.fastq.gz
        
}

export -f run_bbmap

find ${path_input} -type f -name "*R1_001_trimmed_cleaned.fastq.gz" | parallel -j 10 run_bbmap

## track version
 $path_bbmap/reformat.sh --version >> ${path_project_dir}/run_info/tools.txt

# ORIENTING READS
run_vsearch(){
        reads=$1
        sample=$(echo "$reads" | sed 's/_mergedpairs.fastq.gz//')
        sample=$(basename ${sample})
        vsearch --orient ${reads} --db ${orienting_db} --fastqout ${path_output}/oriented/${sample}-oriented.fq 
}

export -f run_vsearch

find ${path_output}/merged/ -type f -name "*mergedpairs.fastq.gz" | parallel -j 10 run_vsearch

## track version
vsearch --version >> ${path_project_dir}/run_info/tools.txt

## FASTQC+MULTIQC for merged and oriented
mkdir ${path_output}/oriented/reports/
find ${path_output}/oriented/ -type f -name "*.fq" | parallel -j 10 "fastqc -o "${path_output}/oriented/reports/" -quiet -t 5"

cd ${path_output}/oriented/reports/
multiqc . --interactive

## track version
fastqc --version >> ${path_project_dir}/run_info/tools.txt
multiqc --version >> ${path_project_dir}/run_info/tools.txt

# The final counts statistics
## Output file
outfile=${path_project_dir}/run_info/read_counts_summary.txt
infile=${path_output}/oriented/reports/multiqc_data/multiqc_general_stats.txt

cd ${path_scripts}/preprocessing
python3 counts_summary.py ${infile} ${outfile} 4

## track version
python --version >> ${path_project_dir}/run_info/tools.txt

# copy the summaries to the run_info file
cp ${path_output}/oriented/reports/multiqc_report.html ${path_project_dir}/run_info/merged_oriented_multiqc_report.html