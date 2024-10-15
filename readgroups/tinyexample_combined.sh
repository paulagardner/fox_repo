#!/bin/bash
#SBATCH --job-name=tinyexample
#SBATCH --output=stdout_tiny.o
#SBATCH --error=stderr_tiny.e
#SBATCH --mail-type=NONE #options: ALL, END, ERROR, etc
#SBATCH --mail-user=paula.gardner@uea.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=3-0 #D-H:M
#SBATCH --partition=compute-64-512


#file where I'm trying to do the Read group thing for only one file, for learning purposes 

module load samtools
module load bwa 

#hard coding to investigate how the pipe anders has written works. 

#JUST print the metadata contained within the config file:
###########################################################
#sed 1d /gpfs/data/bergstrom/foxseq2024/read-group-config.txt \
#| head -1 \
#| while read -r file ID SM LB PL; do \
#echo $file
#done
#this works. Want to get ALL lines from config file? remove the head -1 pipe 

#print out a specified file's read group information  (by its line number in the .txt file). 
#Had to do this as anders was changing file
#############################################################
#sed '1d; 5q; d' /gpfs/data/bergstrom/foxseq2024/read-group-config.txt \
#| while read -r file ID SM LB PL; do \
#echo "$ID $SM $LB $PL" \
#done

#prints:
#227WHLLT4.2.LIS629

#Actually run bwa mem and assign RG data 
sed '1d; 5q; d' /gpfs/data/bergstrom/foxseq2024/read-group-config.txt \
| while read -r file ID SM LB PL; do \
bwa mem -t 8 -T 0 -R "@RG\tID:$ID\tSM:$SM\tLB:$LB\tPL:$PL" \
/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa \
/gpfs/data/bergstrom/foxseq2024/"$file"_1.fq.gz \
/gpfs/data/bergstrom/foxseq2024/"$file"_2.fq.gz \
|samtools view -b -o bamfiles/$ID.bwa.bam ; \
done



#Simply combine what anders did
#################################################
#Soft link read group config file:
#CONFIG=/gpfs/data/bergstrom/foxseq2024/read-group-config.txt
#index=/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa
#Load bwa and samtools:
#module load bwa
#module load samtools
#Loop over files, submit separate run-bwa.sh job for each file - here using the "head" command
#to get only the first line, for testing purposes:
#sed '1d; 5q;d' "$CONFIG" | head -1 | while read file ID SM LB PL; do \
#bwa mem -t 8 -T 0 -R "@RG\tID:$ID\tSM:$SM\tLB:$LB\tPL:$PL" \
#"$index" "$file"_1.fq.gz "$file"_2.fq.gz \
#| samtools view -b -o $ID.bwa.bam ; \
#done
