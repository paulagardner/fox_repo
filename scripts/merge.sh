#!/bin/bash

module load samtools

# Define the input directory
INPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/readgroups/sorted_bamfiles"
OUTPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/data/merged_bamfiles"
LOG_DIR="/gpfs/data/bergstrom/paula/fox_repo/data/merged_bamfiles/slurmout"
TXTPATH="/gpfs/data/bergstrom/foxseq2024/read-group-config.txt"

# Create the logs directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Get the list of unique sample names
SAMPLE_NAMES="$(sed 1d ${TXTPATH} | cut -f 3 | sort | uniq)"
for sample in $SAMPLE_NAMES; do
	echo "Sample: $sample"
	
	# Find all matching files for the current sample
	SAMPLE_FILES=$(ls "${INPUT_DIR}/${sample}"*.bam)
	echo "list of sample files: $SAMPLE_FILES"
	# Check if any files were found for the sample
	if [ -n "$SAMPLE_FILES" ]; then
		# Specify the output file for merged results
		OUTPUT_FILE="${OUTPUT_DIR}/${sample}_merged.bam"
		echo "output file $OUTPUT_FILE"
		# Submit the job with an inline sbatch script
		sbatch -J "merge_${sample}" \
			-o "${LOG_DIR}/${sample}_merge.out" \
			-e "${LOG_DIR}/${sample}_merge.err" \
			--time=2-0 \
			--mem=6G \
			--cpus-per-task=8 \
			--partition=compute-64-512 \
			--wrap="module load samtools && samtools merge -f "$OUTPUT_FILE" $(echo $SAMPLE_FILES)"
		#echo "Submitted merge job for $sample"
		#module add samtools
		##samtools merge -f -o "$OUTPUT_FILE" $(echo $SAMPLE_FILES)
		# Break after the first job submission for testing
		break
	else
		echo "No files found for $sample"
	fi
done

