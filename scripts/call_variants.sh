#!/bin/bash

module load GATK/4.6.0.0

# Define input and output directories
input_dir="/gpfs/data/bergstrom/paula/fox_repo/data/"
mkdir -p "$input_dir/VCF_files/slurmout"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/data/VCF_files"
reference="/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa" #note that this is a different path than other reference genome 
#paths in the rest of the analysis, as we want the index file. check w/ Anders

for bam_file in "$input_dir"/marked_duplicates/*.mdup.bam; do
	# Extract the base name, stripping out both '.merged' and '.bam' extensions if present
	base_name=$(basename "$bam_file" | sed 's/_mdup//; s/\.bam$//')
	echo $base_name
 	echo $reference
	# Define output file paths
 	output_file="${output_dir}/${base_name}.vcf.gz"

		# Check if output file already exists
		if [[ -e "$output_file" ]]; then
    		echo "Output file $output_file already exists. Skipping..."
    		continue
	fi

 	# Submit a Slurm job using --wrap. Use info to make a GVCF, where you use allele specific annottions- good if a site is multi-allelic
 	sbatch --job-name="variant_calling${base_name}" \
        	--error "$output_dir/slurmout/${base_name}.VCF.e" \
        	--output "$output_dir/slurmout/${base_name}.VCF.o" \
         	--time=4-0 \
         	--mem=16G \
         	--cpus-per-task=4 \
         	--partition=compute-64-512 \
         	--wrap="gatk --java-options "-Xmx4g" HaplotypeCaller \
			-R $reference \
			-I $bam_file \
			-O $output_file \
			-ERC GVCF"
	break 
done


