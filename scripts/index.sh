#!/bin/bash

module load samtools/1.16.1

# Define input and output directories
BAMFILES_DIR="/gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates"
mkdir -p "$BAMFILES_DIR/slurmout"

# Loop over each .bam file in the input directory
for bam_file in "$BAMFILES_DIR"/*mdup.bam; do
  # Extract the base name (e.g., "sampleX" from "sampleX.mdup.bam")
  base_name=$(basename "$bam_file" .mdup.bam)

  # Construct the expected .bai file path
  output_index="${bam_file}.bai"

  # Check if the .bai file already exists
  if [[ -e "$output_index" ]]; thden
    echo "Output file $output_index already exists. Skipping..."
    continue
  fi

  # Submit a Slurm job using --wrap
  sbatch --job-name="index_${base_name}" \
    --error "$BAMFILES_DIR/slurmout/${base_name}.indexed.e" \
    --output "$BAMFILES_DIR/slurmout/${base_name}.indexed.o" \
    --time=4-0 \
    --mem=4G \
    --cpus-per-task=4 \
    --partition=compute-64-512 \
    --wrap="samtools index -@ 4 ${bam_file}"
done