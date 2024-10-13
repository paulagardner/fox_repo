#!/bin/bash
#SBATCH --job-name=index_reference
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

bwa index -p ~/bergstromlab/paula/data/mVulVul1.fa /gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa # -p flag lets you easily specify which prefix you want for your files- which, helpfully, means you can send them to whichever directory you'd like
