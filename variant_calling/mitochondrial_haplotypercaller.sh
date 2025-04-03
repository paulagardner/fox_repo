#!/bin/bash

##############NOTE that this should be edited if used again to directly output the files as .gz 
# files, which should be achievable with 

# Set the reference genome and output directory
reference="/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/mtDNA/haplotypecaller"
input_csv="/gpfs/data/bergstrom/paula/fox_repo/all_bam_path_info.csv"  # CSV file with header

# Create output directories if needed
mkdir -p "$output_dir/slurmout"

# Load required modules
module load OpenJDK/jdk-20.0.2

# Read the CSV file and skip the header
{
    read header  # Skip the header line
    while IFS=, read -r sample input_file mean_coverage; do
        # Skip empty lines or if critical fields are missing
        [[ -z "$sample" || -z "$input_file" ]] && continue

        # Use the sample name for the base name
        base_name="${sample}"
        output_file="$output_dir/${base_name}.g.vcf.gz"

        # Check if the output file already exists; if so, skip this file
        if [ -f "$output_file" ]; then
            echo "Output file $output_file already exists. Skipping $input_file."
            continue
        fi

        # Submit the job with sbatch
        sbatch --job-name="haplotypecaller_${base_name}" \
               --error "$output_dir/slurmout/${base_name}.VCF.e" \
               --output "$output_dir/slurmout/${base_name}.VCF.o" \
               --time=6-0 \
               --mem=200G \
               --cpus-per-task=5\
               --partition=compute-64-512 \
               --wrap="java -jar /gpfs/software/ada/gatk/4.6.0.0/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar HaplotypeCaller \
                         -R $reference \
                         -I $input_file \
                         -O $output_file \
                         --sample-ploidy 1 \
                         --intervals chrMT \
                         -ERC GVCF"
    done
} < "$input_csv"


