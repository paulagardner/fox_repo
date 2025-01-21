#!/bin/bash

# Define input and output directories
input_dir="/gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates"
mkdir -p "$input_dir/validation_txt/slurmout"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates/validation_txt"

# Loop over each .bam file in the input directory
for bam_file in "$input_dir"/*.mdup.bam; do
	# Extract the base name, stripping out both '.merged' and '.bam' extensions if present
	base_name=$(basename "$bam_file" | sed 's/_mdup//; s/\.bam$//')

 	# Define output file paths
 	output_file="${output_dir}/ValidateSamfile.${base_name}.txt"

	  # Check if output file already exists
  	if [[ -e "$output_file" ]]; then
    		echo "Output file $output_file already exists. Skipping..."
    		continue
	fi

 	# Submit a Slurm job using --wrap
 	sbatch --job-name="bamfile_validation_${base_name}" \
        	--error "$output_dir/slurmout/${base_name}.validate.e" \
        	--output "$output_dir/slurmout/${base_name}.validate.o" \
         	--time=4-0 \
         	--mem=16G \
         	--cpus-per-task=4 \
         	--partition=compute-64-512 \
         	--wrap="java -jar /gpfs/software/ada/picard/2.24.1/picard.jar ValidateSamFile I=$bam_file O=$output_file M=VERBOSE"
	#break 
done


