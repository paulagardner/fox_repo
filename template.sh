#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=stdout.o
#SBATCH --error=stderr.e
#SBATCH --mail-type=NONE #options: ALL, END, ERROR, etc
#SBATCH --mail-user=paula.gardner@uea.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=360 #D-H:M

echo "test"
