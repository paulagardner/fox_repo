#!/bin/bash

input_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering/autosomes/variant_sites"
reference_fasta="/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa"
mkdir -p "$output_dir/slurmout"

########### Run filtering 

module load bcftools/1.15.1
module load OpenJDK/jdk-20.0.2

for vcf in "$input_dir"/*.vcf.gz; do
  # Make sure it's not a .tbi file
  if [[ "$vcf" != *.tbi ]]; then
    echo "Processing $vcf"
    vcf_name=$(basename "$vcf")
    output="variants_only.$vcf_name"
    sbatch --job-name="filter_${vcf_name}" \
           --error "$output_dir/slurmout/${vcf_name}.e" \
           --output "$output_dir/slurmout/${vcf_name}.o" \
           --time=1-0 \
           --mem=8G \
           --cpus-per-task=7 \
           --partition=compute-64-512 \
           --wrap="zcat $vcf | /gpfs/data/bergstrom/sw/bin/vcf-mask-alleles.pl --mode indels | \
                    /gpfs/data/bergstrom/sw/bin/vcf-mask-alleles.pl --mode third | \
                    bcftools view -a | \
                    bcftools norm -f $reference_fasta | \
                    bcftools view -c 1[:minor] | \
                    bgzip -c - > $output_dir/$output |\
                    bcftools index -f $output_dir/$output"
    #break  # <-- For testing, remove this line when batch processing multiple files
  fi
done
 