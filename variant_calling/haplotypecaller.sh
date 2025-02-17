#!/bin/bash

# Set the reference genome and output directory
reference="/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/out_haplotypecaller"
input_tsv="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/bamfiles.tsv"  # Update this path to your TSV file

# Create output directories if needed
mkdir -p "$output_dir/slurmout"

# Load required modules
# module load GATK/4.6.0.0   # Uncomment if needed
module load OpenJDK/jdk-20.0.2

# Loop over each line (file path) in the TSV file
while IFS= read -r input_file; do
    # Skip empty lines
    [[ -z "$input_file" ]] && continue

    # Get the base name without the .bam extension
    base_name=$(basename "$(dirname "$input_file")")
    # Define the output file name
    output_file="$output_dir/${base_name}.g.vcf"

    # Check if the output file already exists; if so, skip this file
    if [ -f "$output_file" ]; then
        echo "Output file $output_file already exists. Skipping $input_file."
        continue
    fi

    # Submit the job with sbatch
    sbatch --job-name="haplotypecallertest${base_name}" \
           --error "$output_dir/slurmout/${base_name}.VCF.e" \
           --output "$output_dir/slurmout/${base_name}.VCF.o" \
           --time=4-0 \
           --mem=64G \
           --cpus-per-task=1 \
           --partition=compute-64-512 \
           --wrap="java -jar /gpfs/software/ada/gatk/4.6.0.0/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar HaplotypeCaller \
                     -R $reference \
                     -I $input_file \
                     -O $output_file \
                     -ERC GVCF"
done < "$input_tsv"
