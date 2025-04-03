#!/bin/bash

input_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/genotypegvcf"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering"
mkdir -p "$output_dir/slurmout"

####### Make thresholds file in this script so that it's always accurate to what is done 

bam_paths_file="/gpfs/data/bergstrom/paula/fox_repo/all_bam_path_info.csv"
thresholds_file="$output_dir/thresholds-file.txt"

# Create the thresholds-file.txt with header
awk -F',' '
NR > 1 {
    # Extract the sample name and mean coverage
    sample_name = $1
    mean_coverage = $3

    # Calculate the third column (mean_coverage * 1.65)
    third_column_value = mean_coverage * 1.65

    # Print the result in the desired format: sample_name 0 third_column_value
    print sample_name "\t0\t" third_column_value
}' "$bam_paths_file" > "$thresholds_file"



########### Run filtering 

module load bcftools/1.15.1
module load OpenJDK/jdk-20.0.2

for vcf in "$input_dir"/*.vcf.gz; do
  # Make sure it's not a .tbi file
  if [[ "$vcf" != *.tbi ]]; then
    echo "Processing $vcf"
    vcf_name=$(basename "$vcf")
    output="filtered.$vcf_name"
    sbatch --job-name="filter_${vcf_name}" \
           --error "$output_dir/slurmout/${vcf_name}.e" \
           --output "$output_dir/slurmout/${vcf_name}.o" \
           --time=1-0 \
           --mem=4G \
           --cpus-per-task=6 \
           --partition=compute-64-512 \
           --wrap="bcftools +setGT $vcf  -- -t q -n . -i 'FMT/RGQ<20 | GQ<20' | \
                  /gpfs/data/bergstrom/sw/bin/recall-VCF-genotypes.pl --mode sampleDP --DPthresholds $thresholds_file | \
                  bcftools filter -e 'ExcessHet >= 40' -Ou | \
                  bcftools view --include 'F_MISSING < 0.8' -Ou | \
                  bcftools +fill-tags -- -t AN,AC,AF,NS | \
                  bgzip -c - > $output_dir/$output"

    #break  # <-- For testing, remove this line when batch processing multiple files
  fi
done
