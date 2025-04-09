#!/bin/bash

mkdir -p "slurmout"

module load bcftools/1.15.1

sbatch --job-name="subsample" \
        --error "slurmout/downsample.e" \
        --output "slurmout/downsample.o" \
        --time=1-0 \
        --mem=32G \
        --cpus-per-task=5 \
        --partition=compute-64-512 \
        --wrap="bcftools view -c 2[:minor] -Ov Vvulpes_autosomal_variantsonly.vcf.gz | \
                /gpfs/data/bergstrom/sw/bin/sample-lines.pl --leaveHeader 0.1 | \
                bgzip -c - > Vvulpes.c2.samplesites01.vcf.gz | \
                bcftools index -f Vvulpes.c2.samplesites01.vcf.gz"

