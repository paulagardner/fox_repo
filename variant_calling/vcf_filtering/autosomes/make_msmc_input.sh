#!/bin/bash
#SBATCH --output=msmc_files/slurmout/msmc_input_%j.out
#SBATCH --cpus-per-task=3
#SBATCH --mem=2G
#SBATCH --job-name=make_mask_files
#SBATCH --time=1-00:00:00
#SBATCH --array=1-864

#how I made the sample list: ls /gpfs/home/xrq24scu/fox_repo/variant_calling/vcf_filtering/autosomes/het_analysis_filtered/*.chr*.vcf.gz > msmc_input_files.txt



set -euo pipefail

input_list="msmc_input_files.txt"
output_dir="msmc_files/bedfiles"
mkdir -p "$output_dir/slurmout"

module load bcftools/1.22
module load bedtools

# Get file path for this task
file=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" "$input_list")
filename=$(basename "$file")
chr_sample_pair="${filename%.vcf.gz}"
log="$output_dir/${chr_sample_pair}.log"

echo "Processing $file" | tee "$log"

# bcftools commands: Extract CHROM, POS0, and POS fields from VCF using bcftools
# Merge overlapping intervals using bedtools
# Compress the output to .bed.gz using bgzip

bcftools query -f '%CHROM\t%POS0\t%POS\n' "$file" 2>>"$log" | \
    bedtools merge -i - 2>>"$log" | \
    bgzip -c > "$output_dir/${chr_sample_pair}.mask.bed.gz" 2>>"$log"

echo "Finished processing $file" | tee -a "$log"
    
