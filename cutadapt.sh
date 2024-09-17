#!/bin/bash
source activate base
conda init
source ~/.bashrc

# conda environments - change this !!!
#################################################################
conda_envs_dir=/home/povp/conda_envs
r_env=/home/povp/conda_evs/R_v4_env
###################################################################

###############################################################
# personalize these paths!!!!
path="/home/povp/16S/knihovna9/"
path_output="/home/povp/spracovane_data/knihovna9/"
path_base="/home/povp/"
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

# vytvor zlozky
mkdir ${path_output}kvalita_raw
mkdir ${path_output}trimmed
mkdir ${path_output}kvalita_trimmed
cd ${path}

# kvalita raw readov
#cd $conda_envs_dir
conda activate /home/povp/conda_envs/quality

fastqc *.fastq.gz -o ${path_output}kvalita_raw --quiet

cd ${path_output}kvalita_raw
multiqc . --interactive
cd multiqc_data

# statistika
awk -F'\t' 'NR==1 {print "SampleName\tRawReads"} {print $1"\t"$NF}' multiqc_general_stats.txt > ${path_output}raw_reads.csv

# adapters + primers trimming 
cd $path

for read1 in *R1_001.fastq*
do
	read2=$(echo $read1| sed 's/R1_/R2_/')
	command1="cutadapt "
	cd $path_base
	for item in $(cat microbiome_primers.csv) 
	do
		FWD=$(echo $item | cut -f2 -d ",")
		REV=$(echo $item | cut -f5 -d ",")

		command1+="-g ^${FWD} -G ^${REV} "
	done
	command1+="-a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG --discard-untrimmed -j 10 -o ${path_output}/trimmed/t2_${read1} -p ${path_output}/trimmed/t2_${read2} \
	${path}/${read1} ${path}/${read2}"

	eval $command1
	
done

# kvalita trimmed readov

cd $path_output/trimmed
fastqc *.fastq.gz -o ${path_output}kvalita_trimmed --quiet

cd ${path_output}kvalita_trimmed
multiqc . --interactive

# statistika
cd multiqc_data
awk -F'\t' 'NR==1 {print "TrimmedReads"} {print $NF}' multiqc_general_stats.txt > ${path_output}trimmed_reads.csv

cd ${path_output}
paste -d '\t' raw_reads.csv trimmed_reads.csv | awk -F'\t' 'BEGIN {OFS="\t"} {print $0}' > output.csv
sed '2d' output.csv > statistics.csv

rm raw_reads.csv
rm trimmed_reads.csv
rm output.csv











