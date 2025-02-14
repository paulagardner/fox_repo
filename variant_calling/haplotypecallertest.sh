#!/bin/bash

reference="/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa" 
input_file="/gpfs/data/bergstrom/paula/fox_repo/Hasselgren_genomes/eager_results/deduplication/15075/15075_rmdup.bam"
output_file="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/testfile.g.vcf.gz"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling"
base_name=$(basename "$input_file" .bam)
mkdir -p "$output_dir/slurmout"

#module load GATK/4.6.0.0

module load OpenJDK/jdk-20.0.2


sbatch --job-name="haplotypecallertest${base_name}" \
        --error "$output_dir/slurmout/${base_name}.VCF.e" \
        --output "$output_dir/slurmout/${base_name}.VCF.o" \
        --time=4-0 \
        --mem=64G  \
        --cpus-per-task=1 \
        --partition=compute-64-512 \
        --wrap="java -jar /gpfs/software/ada/gatk/4.6.0.0/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar HaplotypeCaller \
		-R $reference \
		-I $input_file \
		-O $output_dir/testfile.g.vcf.gz \
		-ERC GVCF"