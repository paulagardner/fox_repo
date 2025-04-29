#!/bin/bash
#SBATCH --job-name=vcf_processing  # Job name
#SBATCH --output=slurm-%j.out      # Standard output log
#SBATCH --error=slurm-%j.err       # Standard error log
#SBATCH --time=02:00:00            # Walltime
#SBATCH --mem=8G                   # Memory allocation
#SBATCH --cpus-per-task=4          # Number of CPUs per task
INPUT_VCF_DIR="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering/autosomes"
OUTPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/MSMC/splitbysample"
mkdir -p "${OUTPUT_DIR}/logs"

#module load bcftools/1.15.1

#for VCF in "$INPUT_VCF_DIR"/*.vcf.gz; do

#    BASENAME=$(basename "$VCF" .vcf.gz)
#    for SAMPLE in $(bcftools query -l "$VCF"); do
#        SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
#        mkdir -p "$SAMPLE_DIR"

#        OUTFILE="${SAMPLE_DIR}/${SAMPLE}_${BASENAME}.vcf.gz"

#        sbatch --job-name="vcf_${SAMPLE}" \
#               --output="${OUTPUT_DIR}/logs/${SAMPLE}_${BASENAME}_%j.out" \
#               --error="${OUTPUT_DIR}/logs/${SAMPLE}_${BASENAME}_%j.err" \
#               --time=01:00:00 \
#               --mem=4G \
#               --cpus-per-task=1 \
#               --wrap="bcftools view -c 1 -s $SAMPLE -Oz -o $OUTFILE $VCF && \
#                       bcftools index $OUTFILE"
#    done
#    break
#done


module load python/3.10.0

# Loop through all .vcf.gz files in subdirectories of $OUTPUT_DIR
for file in "$OUTPUT_DIR"/*/*.gz; do
  # Extract the base filename without the .gz extension
  base=$(basename "$file" .gz)

  # Extract the name of the subdirectory containing the file
  subdir=$(basename "$(dirname "$file")")

  # Define the output directory where results will be saved
  outdir="$OUTPUT_DIR/newdir"  # Assuming output goes to $OUTPUT_DIR/newdir
  mkdir -p "$outdir/$subdir"  # Create subdirectory within newdir if it doesn't exist

  # Extract the chromosome name from the first non-header line in the VCF file
  chrom=$(zcat "$file" | grep -v '^#' | head -n1 | cut -f1)

  # Define the output filenames for the mask (BED format) and filtered VCF
  mask_out="$outdir/$subdir/${base}.mask.bed.gz"
  vcf_out="$outdir/$subdir/${base}.filtered.vcf"

  # Now, wrap only the Python script part to be run by srun --wrap
  srun --wrap="python3 /gpfs/data/bergstrom/paula/fox_repo/MSMC/vcfAllSiteParser.py \"$chrom\" \"$mask_out\" \"$file\" > \"$vcf_out\""
done
