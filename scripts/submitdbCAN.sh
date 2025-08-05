#!/bin/bash
#
#$ -q all.q
#$ -cwd
#$ -pe smp 20

module load dbCAN/5.1.2

for file in *.faa; do
	base=${file/.faa}
	if [ ! -d ${base}_CAZymes ]; then
		run_dbcan CAZyme_annotation --input_raw_data ${file} --output_dir ${base}_CAZymes --db_dir /Storage/databases/dbCAN_v5.1.2 --mode protein --threads $NSLOTS
	fi
done
