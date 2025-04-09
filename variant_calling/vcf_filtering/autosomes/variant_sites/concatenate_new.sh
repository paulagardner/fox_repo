#!/bin/bash

mkdir -p "slurmout"

module load bcftools/1.15.1

sbatch --job-name="concatenate" \
        --error "slurmout/concat.e" \
        --output "slurmout/concat.o" \
        --time=1-0 \
        --mem=32G \
        --cpus-per-task=4 \
        --partition=compute-64-512 \
        --wrap="ls -1 *chr*.vcf.gz | sort -V > files-to-concat.txt | \
                bcftools concat -f files-to-concat.txt -Oz -o 
                /gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering/autosomes/variant_sites/autosomal_variantsonly.vcf.gz | \
                bcftools index -f \
                /gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering/autosomes/variant_sites/autosomal_variantsonly.vcf.gz"

