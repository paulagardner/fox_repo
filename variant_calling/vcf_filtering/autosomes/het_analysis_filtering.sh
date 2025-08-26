#!/bin/bash
#SBATCH --job-name=mask_vcf
#SBATCH --output=het_analysis_filtered/logs/mask_chr_%A_%a.out
#SBATCH --error=het_analysis_filtered/logs/mask_chr_%A_%a.err
#SBATCH --array=0-15       # One job per chromosome file
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00

# Optional: load modules or activate conda environment
# module load bcftools
# source activate your_env
module load bcftools

mkdir -p het_analysis_filtered/logs


# List all VCFs matching the chromosome pattern
vcf_files=(filtered.out.chr*.vcf.gz)
vcf_file="${vcf_files[$SLURM_ARRAY_TASK_ID]}"
chrom=$(basename "$vcf_file" .vcf.gz | cut -d. -f3)

echo "[$(date)] Starting processing for chromosome: $chrom"
echo "VCF file: $vcf_file"
echo "Using $SLURM_CPUS_PER_TASK CPU threads"



#bcftools query -l : will list the samples from the samples fields, which we use to loop
#then, we pull up the whole file with bcftools view
# remove insertions + deletions with anders' custom script
# re-format with bcftools norm: "left align and normalize indels"
# bcftools view --exclude-uncalled : exclude sites without a called genotype
# bcftools annotate -x INFO,QUAL,FILTER : remove info from these fields from the output.
###we keep the fields in to still have the file be the correct format
#bcftools annotate -x ^FMT/GT .... etc: removes all fields from the VCF except the genotype (GT)
### field in the FORMAT section, then compresses the result and writes it to a file. -z compresses output



# Loop over each sample in the VCF
for sample in $(bcftools query -l "$vcf_file"); do
    echo "→ Processing sample: $sample"

    out="het_analysis_filtered/${chrom}.${sample}.masked.vcf.gz" 

    bcftools view --threads $SLURM_CPUS_PER_TASK -s "$sample" "$vcf_file" \

    | /gpfs/data/bergstrom/sw/bin/vcf-mask-alleles.pl --mode indels \
    | bcftools view -a -Ou --threads $SLURM_CPUS_PER_TASK \
    | bcftools norm -f /gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa -Ou --threads $SLURM_CPUS_PER_TASK \
    | bcftools view --exclude-uncalled -Ou --threads $SLURM_CPUS_PER_TASK \
    | bcftools annotate -x INFO,QUAL,FILTER -Ou --threads $SLURM_CPUS_PER_TASK \
    | bcftools annotate -x ^FMT/GT -Oz --threads $SLURM_CPUS_PER_TASK -o "$out"

    echo "✓ Done with $sample"
done

echo "[$(date)] Finished processing chromosome: $chrom"