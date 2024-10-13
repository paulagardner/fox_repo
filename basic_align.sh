#!/bin/bash
#SBATCH --job-name=basic_bwa
#SBATCH --output=stdout.o
#SBATCH --error=stderr.e
#SBATCH --mail-type=NONE #options: ALL, END, ERROR, etc
#SBATCH --mail-user=paula.gardner@uea.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=360 #D-H:M

module load bwa/0.7.17
module load samtools/1.16.1

bwa mem /gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa \
/gpfs/data/bergstrom/foxseq2024/S080738/S080738_EKDN240031781-1A_22MGFYLT3_L2_1.fq.gz \
/gpfs/data/bergstrom/foxseq2024/S080738/S080738_EKDN240031781-1A_22MGFYLT3_L2_2.fq.gz \
> /gpfs/data/bergstrom/paula/data/test-output.sam
