#!/bin/bash

TXTPATH="/gpfs/data/bergstrom/foxseq2024/read-group-config.txt"
#Load bwa and samtools:
module load samtools

mkdir -p sorted_bamfiles
mkdir -p sorted_bamfiles/slurmout
mkdir -p sorted_bamfiles/temp

for file in ~/fox_repo/readgroups/readgroup_bamfiles/*.bam; do
    # Extract the prefix (basename without extension)
    prefix=$(basename "$file" .bam) #putting .bam at the end EXCLUDES it from prefix, which will allow you to put sorted BEHIND a .bam extension
    echo "$prefix"
    # Submit the job using sbatch
    sbatch -J sort_samtools --time 3-0 --mem=12G --cpus-per-task=4 --partition=compute-64-512 \
    -o sorted_bamfiles/slurmout/"$prefix".sorted.o -e sorted_bamfiles/slurmout/"$prefix".sorted.e \
    --wrap="samtools sort -@ 4 -m 2G -T sorted_bamfiles/temp/'$prefix'.temp --reference /gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa '$file' -o sorted_bamfiles/'$prefix'.sorted.bam" 
    # include a break after processing the first file for testing purposes
    #break
done





