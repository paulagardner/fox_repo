#!/bin/bash

mkdir -p "slurmout"

module load bcftools/1.15.1

sbatch --job-name="subsample" \
        --error "slurmout/redfoxsubsample.e" \
        --output "slurmoutredfoxsubsample.o" \
        --time=1-0 \
        --mem=32G \
        --cpus-per-task=5 \
        --partition=compute-64-512 \
        --wrap="bcftools view -S samples_to_keep.txt -Ou autosomal_variantsonly.vcf.gz | \
                bcftools view -c 1[:minor] -Oz > Vvulpes_autosomal_variantsonly.vcf.gz | \
                bcftools index Vvulpes_autosomal_variantsonly.vcf.gz"

