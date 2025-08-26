#!/bin/bash

input=$1    # The input VCF
sample=$2   # The sample ID to extract
output=$3   # The output VCF to be created

# Process the VCF for the given sample
bcftools view -s $sample $input \
| /gpfs/data/bergstrom/sw/bin/vcf-mask-alleles.pl --mode indels \
| bcftools view -a -Ou \
| bcftools norm -f /gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa -Ou \
| bcftools view --exclude-uncalled -Ou \
| bcftools annotate -x INFO,QUAL,FILTER -Ou \
| bcftools annotate -x ^FMT/GT -Oz -o $output

#bcftools query -l : will list the samples from the samples fields, which we use to loop
#then, we pull up the whole file with bcftools view
# remove insertions + deletions with anders' custom script
# re-format with bcftools norm: "left align and normalize indels"
# bcftools view --exclude-uncalled : exclude sites without a called genotype
# bcftools annotate -x INFO,QUAL,FILTER : remove info from these fields from the output.
###we keep the fields in to still have the file be the correct format
#bcftools annotate -x ^FMT/GT .... etc: removes all fields from the VCF except the genotype (GT)
### field in the FORMAT section, then compresses the result and writes it to a file. -z compresses output
