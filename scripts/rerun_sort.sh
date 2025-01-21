#!/bin/bash

INPUT_FILE="/gpfs/data/bergstrom/paula/fox_repo/data/rerun_sort.txt"
INPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/readgroups/readgroup_bamfiles/"
OUTPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/data/sorted_bamfiles"
REFERENCE_GENOME="/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa"
TEMP_DIR="/gpfs/scratch/xrq24scu"
module load samtools

# Ensure output directories exist
mkdir -p "$OUTPUT_DIR/slurmout"

# Loop through each line in the input file without using IFS= read -r
for line in $(cat "$INPUT_FILE"); do
    # Ensure the line is treated as a full file path (to avoid issues with relative paths)
    full_path="${INPUT_DIR}${line}"

    # Extract the prefix (basename without extension) for each file in the list
    prefix=$(basename "$line" .bam)

    # Submit the job using sbatch
    sbatch -J sort_samtools --time=4-0 --mem=20G --cpus-per-task=4 --partition=compute-64-512 \
    -o "$OUTPUT_DIR/slurmout/$prefix.sorted.o" \
    -e "$OUTPUT_DIR/slurmout/$prefix.sorted.e" \
    --wrap="samtools sort -@ 6 -m 2G -T $TEMP_DIR/$prefix.temp $full_path -o $OUTPUT_DIR/$prefix.sorted.bam"

    # Uncomment the following line to process only the first file for testing
    # break
done

