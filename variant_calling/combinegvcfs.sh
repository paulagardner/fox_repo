#!/bin/bash

# Set the reference genome and output directory
reference="/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/combinegvcfs"
gvcf_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/adjustednames_haplotypecaller"  # Directory containing GVCF files

# Create output directories if needed
mkdir -p "$output_dir/slurmout"

# Load required modules
module load OpenJDK/jdk-20.0.2

# Construct the --variant arguments, ensuring only .g.vcf files are selected
# and NOT their indexes 
VARIANT_ARGS=""
for file in "$gvcf_dir"/*.g.vcf; do
    [[ -f "$file" ]] && VARIANT_ARGS+="--variant $file "
done


### -Xmx argument comes from wrapper file instructions: https://gatk.broadinstitute.org/hc/en-us/articles/360035531892-GATK4-command-line-syntax
## and troubleshooting on memory: https://gatk.broadinstitute.org/hc/en-us/community/posts/360074224671-GenomicsDBImport-running-out-of-memory

# Submit the job with sbatch
sbatch --job-name="combinegvcfs" \
        --error "$output_dir/slurmout/allchromcohort.VCF.e" \
        --output "$output_dir/slurmout/allchromcohort.VCF.o" \
        --time=3-0 \
        --mem=200G \
        --cpus-per-task=4 \
        --partition=compute-64-512 \
        --wrap="java -jar -Xmx190g -XX:ConcGCThreads=4  /gpfs/software/ada/gatk/4.6.0.0/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar CombineGVCFs \
                --tmp-dir /gpfs/scratch/xrq24scu \
                -R $reference \
                $VARIANT_ARGS \
                -O $output_dir/allchromcohort.g.vcf.gz"

