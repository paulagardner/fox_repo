#!/bin/bash

INPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/readgroups/readgroup_bamfiles"
OUTPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/data/sorted_bamfiles"
REFERENCE_GENOME="/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa"
TEMP_DIR="/gpfs/scratch/xrq24scu"

module load samtools

# Ensure output directories exist
mkdir -p "$OUTPUT_DIR/slurmout"

for file in $INPUT_DIR/*.bam; do
    # Extract the prefix (basename without extension)
    prefix=$(basename "$file" .bam)

    # Submit the job using sbatch
    sbatch -J sort_samtools --time 4-0 --mem=16G --cpus-per-task=4 --partition=compute-64-512 \
    -o "$OUTPUT_DIR/slurmout/$prefix.sorted.o" -e "$OUTPUT_DIR/slurmout/$prefix.sorted.e" \
    --wrap="samtools sort -@ 4 -m 2G -T $TEMP_DIR/$prefix.temp --reference $REFERENCE_GENOME $file -o $OUTPUT_DIR/$prefix.sorted.bam"
    # Uncomment the following line to process only the first file for testing
    #break
done
