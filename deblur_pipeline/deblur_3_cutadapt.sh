source /opt/anaconda3/etc/profile.d/conda.sh

conda activate /home/povp/conda_envs/qiime2-amplicon-2024.5
export TMPDIR=/home/povp/data

path=/home/povp/16S/
path_vysledky=/home/povp/16S_deblur/
path_base=/home/povp/
path_bbmap=/home/povp/bbmap/

#CUTADAPT
for lib_name in knihovna2 knihovna3 knihovna4 knihovna5 knihovna6 knihovna7 knihovna8 knihovna9
	do	
	mkdir ${path_vysledky}/${lib_name}/trimmed
        cd ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk
        for read1 in *R1_001.fastq*
	do
	read2=$(echo $read1| sed 's/R1_/R2_/')
	sample1=$(echo $read1| sed 's/noadaptors_nophix_//')
	sample2=$(echo $read2| sed 's/noadaptors_nophix_//')
	command1="cutadapt "
	cd $path_base
	for item in $(cat microbiome_primers.csv)
	do
	FWD=$(echo $item | cut -f2 -d ",")
	REV=$(echo $item | cut -f5 -d ",")
        command1+="-g ^${FWD} -G ^${REV} "
	done
	command1+="--discard-untrimmed -j 5 -o ${path_vysledky}/${lib_name}/trimmed/trimmed_${sample1} -p ${path_vysledky}/${lib_name}/trimmed/trimmed_${sample2} \
	${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/${read1} ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/${read2}"
	eval $command1
done
done

