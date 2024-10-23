#!/bin/bash

#read command line arguments
file=$1		#prefix for paired FASTQ files, the below command assumes names are followed by "_1.fq.gz" and "_2.fq.gz"
ID=$2
SM=$3
LB=$4
PL=$5

#hard code bwa index location
index=/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa
output=$(basename $file)

bamfilesdir=readgroup_bamfiles
mkdir -p "$bamfilesdir"

#map reads with bwa, and convert output to BAM on the fly using samtools
bwa mem -t 8 -T 0 -R "@RG\tID:$ID\tSM:$SM\tLB:$LB\tPL:$PL" $index "$file"_1.fq.gz "$file"_2.fq.gz | samtools view -b -o "$bamfilesdir"/$output.bwa.bam
