#!/bin/bash
eval "$(conda shell.bash hook)"

# variables for conda environment
amplicon_16s_env=$1

# paths for data storage
path_input=$2
path_project_dir=$3
path_output=${path_project_dir}/preprocessed

# activate environment
conda activate ${amplicon_16s_env}

# info for tools and versions txt
echo -e "\ncutadapt.sh:" >> ${path_project_dir}/run_info/tools.txt

# create output dirs
mkdir ${path_output}
mkdir ${path_output}/reports/

# adapters + primers trimming 

for read1 in ${path_input}/*R1_001.fastq.gz
do
	read2=$(echo $read1| sed 's/R1_/R2_/')
    sample_name=$(basename "${read1}" | sed 's/_R1_001.fastq.gz//')
	command1="cutadapt "
	for item in $(cat microbiome_primers.csv) 
	do
		FWD=$(echo $item | cut -f2 -d ",")
		REV=$(echo $item | cut -f5 -d ",")

		command1+="-g ^${FWD} -G ^${REV} "
	done
	command1+="-a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
    -a 'A{10}' -a 'G{10}' -g 'A{10}' -g 'G{10}' \
    --discard-untrimmed -j 10 \
    -o ${path_output}/${sample_name}_R1_001_trimmed.fastq.gz \
    -p ${path_output}/${sample_name}_R2_001_trimmed.fastq.gz \
	${read1} ${read2}"

	eval $command1
	
done

## track version
echo "cutadapt " $(cutadapt --version) >> ${path_project_dir}/run_info/tools.txt

# FASTQC+MULTIQC for trimmed report
find ${path_output} -type f -name "*.fastq.gz" | parallel -j 10 "fastqc -o "${path_output}/reports/" -quiet -t 5"
cd ${path_output}/reports/
multiqc . --interactive

## track version
fastqc --version >> ${path_project_dir}/run_info/tools.txt
multiqc --version >> ${path_project_dir}/run_info/tools.txt

# The final counts statistics
## Output file
outfile=${path_project_dir}/run_info/read_counts_summary.txt
infile=${path_output}/reports/multiqc_data/multiqc_general_stats.txt

cd ${path_scripts}/preprocessing
python3 counts_summary.py ${infile} ${outfile} 2

## track version
python --version >> ${path_project_dir}/run_info/tools.txt

# copy the summaries to the run_info file
cp ${path_output}/reports/multiqc_report.html ${path_project_dir}/run_info/trimmed_multiqc_report.html
