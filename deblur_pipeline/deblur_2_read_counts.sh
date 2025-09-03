#!/bin/bash

source /opt/anaconda3/etc/profile.d/conda.sh

conda activate /home/povp/conda_envs/qiime2-amplicon-2024.5
export TMPDIR=/home/povp/data

path=/home/povp/16S/
path_vysledky=/home/povp/16S_deblur/
path_base=/home/povp/
path_bbmap=/home/povp/bbmap/


cd $path

for lib_name in knihovna2 knihovna3 knihovna4 knihovna5 knihovna6 knihovna7 knihovna8 knihovna9
do
	cd ${path}/${lib_name}/
	echo -e "sampleid\tbases_in_R1\treads_in_R1\tavg_read_length_in_R1" > ${path_vysledky}/${lib_name}/readcount_raw_${lib_name}.tsv
	for i in $(ls *_R1_* | cut -f 1)
	do
		echo -en "$i\t" >> ${path_vysledky}/${lib_name}/readcount_raw_${lib_name}.tsv; \
			zcat ${i} | paste - - - - | cut -f 2 | perl -ne 'END {if ($. == 0) {print "0\t0\t0\n"} else {print "$c\t$.\t", sprintf("%.0f", $c/$.), "\n"}} $c+=(length($_)-1)' - >> ${path_vysledky}/${lib_name}/readcount_raw_${lib_name}.tsv; \
		done
	done

