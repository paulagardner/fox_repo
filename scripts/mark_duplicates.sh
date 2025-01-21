#!/bin/bash

# Define input and output directories
input_dir="/gpfs/data/bergstrom/paula/fox_repo/data/merged_bamfiles"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates"
mkdir -p "$output_dir/slurmout"

# Loop over each .bam file in the input directory
for bam_file in "$input_dir"/*.bam; do
	# Extract the base name, stripping out both '.merged' and '.bam' extensions if present
	base_name=$(basename "$bam_file" | sed 's/_merged//; s/\.bam$//')

 	# Define output file paths
	output_bam="${output_dir}/${base_name}.mdup.bam"
 	metrics_file="${output_dir}/metrics-MarkDuplicates.${base_name}.txt"

	  # Check if output file already exists
  	if [[ -e "$output_bam" ]]; then
    		echo "Output file $output_bam already exists. Skipping..."
    		continue
	fi

 	# Submit a Slurm job using --wrap
 	sbatch --job-name="mark_duplicates_${base_name}" \
        	--error "$output_dir/slurmout/${base_name}.mdup.e" \
        	--output "$output_dir/slurmout/${base_name}.mdup.o" \
         	--time=4-0 \
         	--mem=16G \
         	--cpus-per-task=4 \
         	--partition=compute-64-512 \
         	--wrap="java -jar /gpfs/software/ada/picard/2.24.1/picard.jar MarkDuplicates I=$bam_file O=$output_bam M=$metrics_file"
done


