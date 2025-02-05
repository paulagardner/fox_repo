import glob

#module load python/anaconda/2024.06/3.12.4

### this matching pattern assumes that eager_results is the name of the results dir.
### Could rewrite this to reference what is assigned in the wrappper, but currently using
### one wrapper per directory, so this depends on all the wrappers having matching
### results names.


from pathlib import Path
import csv

# Define the base directory
base_dir = Path("/gpfs/data/bergstrom/paula/fox_repo")

# Define the relative glob pattern
pattern = "**/eager_results/deduplication/**/*_rmdup.bam"

# Use the .glob() method to find matching files
files = base_dir.glob(pattern)

# Specify the output CSV file path
output_csv = "all_bam_path_info.csv"

# Open the CSV file for writing
with open(output_csv, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write the header row with 'sample' first, then 'path'
    writer.writerow(['sample_ID', 'path'])
    
    # For each matching file, extract the sample and write the row.
    for file in files:
        # Extract the sample from the filename by taking the substring before the first underscore.
        # sample = file.name.split('_')[0]
        sample = file.parent.name
        writer.writerow([sample, str(file)])
