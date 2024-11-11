#!/bin/bash

# Load the picard module
#module load picard/2.24.1

#module show picard #to find the path
"java -jar /gpfs/software/ada/picard/2.24.1/picard.jar"


# Define input and output directories
input_dir="/gpfs/data/bergstrom/paula/fox_repo/data/merged_bamfiles"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates"
mkdir -p "$output_dir/slurmout"
# Loop over each .bam file in the input directory
for bam_file in "$input_dir"/*.bam; do
  # Extract the base name (e.g., "sampleX" from "sampleX.merged.bam")
	base_name=$(basename "$bam_file" .merged.bam)
  # Define output file paths
	output_bam="${output_dir}/${base_name}.mdup.bam"
	metrics_file="${output_dir}/metrics-MarkDuplicates.${base_name}.txt"

  # Submit a Slurm job using --wrap
	sbatch --job-name="mark duplicates" \
		--error "$output_dir/slurmout/${base_name}.e" \
		--output "$output_dir/slurmout/${base_name}.o" \
		--time=4-0 \
		--mem=4G \
		--cpus-per-task=4 \
		--partition=compute-64-512 \
		--wrap="java -jar /gpfs/software/ada/picard/2.24.1/picard.jar MarkDuplicates I=$bam_file O=$output_bam M=$metrics_file"
done
