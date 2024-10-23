#!/bin/bash


#path to read groups config file:
TXTPATH="/gpfs/data/bergstrom/foxseq2024/read-group-config.txt"

#Load bwa and samtools:
module load bwa
module load samtools

#Loop over files, submit separate run-bwa.sh job for each file - here using the "head" command to get only the first line, for testing purposes:
slurmdir=slurm_output
mkdir -p "$slurmdir" #https://stackoverflow.com/questions/2743673/mkdir-error-in-bash-script/2743818#2743818


#removing |head -1 | pipe from after sed command, as I want to run these all!
sed 1d "$TXTPATH" | while read file ID SM LB PL; do output=$(basename "$file") ; \
sbatch -J batch_align --time=3-0 --mem=10G --cpus-per-task=9 --partition=compute-64-512 \
-o "$slurmdir"/out_bwa."$output".o -e "$slurmdir"/out_bwa."$output".e \
--wrap="./run-bwa.sh /gpfs/data/bergstrom/foxseq2024/$file $ID $SM $LB $PL"; \
done

#each bwa command needs 10GB, but the runtime will depend on the size of the input file. Asking for three days might be reasonable (--time=3-0)
