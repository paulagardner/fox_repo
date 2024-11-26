#!/bin/bash

#columns we're trying to make, in order:
#genetic sex, total reads, fraction mapped, 
#coverage, percent duplicates
#not sure whether to include FASTQs and .BAMs? 

###########Anders' script for number of reads and fraction mapped 
echo "##########NUMBER OF READS AND FRACTION MAPPED##########"
cat /gpfs/data/bergstrom/foxseq2024/samples-gdoc-order.txt \
| while read sample; do (echo $sample; grep "primary mapped" /gpfs/data/bergstrom/paula/fox_repo/data/reads_validation/$sample.flagstat \
| awk '{ print $1"\t"$6}' \
| tr -d "(" ) \
| paste - -; done #> summary.flagstat.txt #leave this commented to simply print to 

#| while read sample; do echo $sample; done  #this line prints all the names in the file


################Percent duplicates:
echo "##########PERCENT DUPLICATES##########" 
cat /gpfs/data/bergstrom/foxseq2024/samples-gdoc-order.txt | while read sample; do cat /gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates/metrics-MarkDuplicates.$sample.txt \
| grep -A1 "^LIBRARY" \
| cut -f 1,9 \
| sed 1d; done


##################Coverage: 
echo "##########AVERAGE COVERAGE##########" 
cat  /gpfs/data/bergstrom/foxseq2024/samples-gdoc-order.txt | while read sample; do cat /gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates/validation_txt/WgsMetrics.$sample.mdup.txt \
| grep -A1 "^GENOME_TERRITORY" \
| cut -f 2 \
| sed 1d \
| awk -v sample="$sample" '{print sample, $0}'; done
