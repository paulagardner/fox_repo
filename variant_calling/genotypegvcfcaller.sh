#!/bin/bash

# Set the reference genome and output directory
reference="/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/genotypegvcf"
input_file="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/combinegvcfs/allchromcohort.g.vcf.gz"

# Create output directories if needed
mkdir -p "$output_dir/slurmout"

# Load required modules
module load OpenJDK/jdk-20.0.2

# Submit the job with sbatch
sbatch --job-name="gvcfcaller" \
        --error "$output_dir/slurmout/out.chr%a.vcf.e" \
        --output "$output_dir/slurmout/out.chr%a.vcf.o" \
        --time=2-0 \
        --mem=25G \
        --cpus-per-task=4 \
        --partition=compute-64-512 \
        --array=1-16 \
        --wrap="java -jar /gpfs/software/ada/gatk/4.6.0.0/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar GenotypeGVCFs \
                --include-non-variant-sites \
                -L chr\${SLURM_ARRAY_TASK_ID} \
                -R $reference \
                -V $input_file \
                --output $output_dir/out.chr\${SLURM_ARRAY_TASK_ID}.vcf.gz"

#However, inside the --wrap command, you must use ${SLURM_ARRAY_TASK_ID} because it runs within the job environment.
