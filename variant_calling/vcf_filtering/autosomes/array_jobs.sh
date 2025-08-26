#!/bin/bash

mkdir -p het_analysis_filtered/logs

# Dynamically determine the chromosome numbers from the files
chromosomes=$(ls filtered.out.chr*.vcf.gz | while read vcf_file; do
    chrom=$(basename "$vcf_file" .vcf.gz | cut -d. -f3)
    echo "$chrom"
done | sort -n | uniq)

# Loop over all chromosomes and samples
for chrom in $chromosomes; do
    bcftools query -l filtered.out.chr${chrom}.vcf.gz | while read sample; do
        sbatch --mem=1G --time=360 \
        -o het_analysis_filtered/logs/out_${sample}.chr${chrom}.o \
        -e het_analysis_filtered/logs/out_${sample}.chr${chrom}.e \
        --wrap="./make-single-sample-VCF.sh filtered.out.chr${chrom}.vcf.gz $sample ${sample}.chr${chrom}.vcf.gz"
    done
done