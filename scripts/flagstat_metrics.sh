#!/bin/bash
#### script to 
module load samtools 

READPAIR_FILE="/gpfs/data/bergstrom/foxseq2024/count-reads/read-pairs-per-sample-foxseq2024.txt"
BAMFILES_DIR="/gpfs/data/bergstrom/paula/fox_repo/data/merged_bamfiles"
FLAGSTAT_DIR="/gpfs/data/bergstrom/paula/fox_repo/data/reads_validation"
mkdir -p "$FLAGSTAT_DIR/slurmout"
#have awk read in the read pairs file BUT do not get the header.
awk 'NR >1'  "$READPAIR_FILE" | while read -r line; do
	#assign variables to the first and second columns, representing the sample and 
	#number of reads
	sample="$(echo "$line" | awk '{print $1}')"
	orig_readpairs="$(echo "$line" | awk '{print $2}')"
	matched_file=$(ls "$BAMFILES_DIR" | grep "^$sample" | head -n 1)
	#check for matches between the files in the merged_bamfiles directory, and the readpair file
	if [ -n "$matched_file" ]; then  
		echo "found a match for $sample, with $orig_readpairs read pairs"
		echo "$matched_file" 
		sbatch	--job-name="$sample_flagstat" \
			--out="$FLAGSTAT_DIR/slurmout/$sample.o" \
			--error="$FLAGSTAT_DIR/slurmout/$sample.e" \
			--time=02:00:00 \
			--mem=4G \
			--cpus-per-task=4 \
			--partition=compute-64-512 \
			--wrap="echo 'original readpairs: $orig_readpairs' > $FLAGSTAT_DIR/$sample.flagstat && \
				samtools flagstat $BAMFILES_DIR/$matched_file >> $FLAGSTAT_DIR/$sample.flagstat"
	else 
		echo "failed to match $sample"
	fi
	#break
done
