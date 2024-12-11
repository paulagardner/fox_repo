#!/bin/bash

# Path to read groups config file
TXTPATH="/gpfs/data/bergstrom/foxseq2024/read-group-config.txt"

# Load bwa and samtools
module load bwa
module load samtools

# Directory to store SLURM output
slurmdir=slurm_output
mkdir -p "$slurmdir"

# Hard-code BWA index location
index=/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa
bamfilesdir=readgroup_bamfiles
mkdir -p "$bamfilesdir"

# Read the config file (skipping the header) and process only the first line
sed 1d "$TXTPATH" | while read file ID SM LB PL; do 
  output=$(basename "$file")

  # Submit the job to SLURM
  sbatch -J batch_align \
    --time=3-0 --mem=10G --cpus-per-task=9 --partition=compute-64-512 \
    -o "$slurmdir"/out_bwa."$output".o -e "$slurmdir"/out_bwa."$output".e \
    --wrap="bwa mem -t 8 -T 0 -R '@RG\tID:$ID\tSM:$SM\tLB:$LB\tPL:$PL' $index \
    /gpfs/data/bergstrom/foxseq2024/${file}_1.fq.gz /gpfs/data/bergstrom/foxseq2024/${file}_2.fq.gz | \
    samtools view -b -o $bamfilesdir/$output.bwa.bam"

  # Exit after processing the first file
  break
done
