#!/bin/bash

# Files
ENA_FILE="filereport_read_run_PRJNA1158265_tsv.txt"
SAMPLES_FILE="samples_of_interest.txt"
OUTPUT_DIR="downloaded_files"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"


# Define the columns to use for this run (you can modify this as needed)
ENA_columns = {
    'ENA_file_id': 'library_name',  # Column in the ENA file for library name (generalized)
    'ftp_links': 'fastq_ftp'        # Column in the ENA file for the fastq URLs (generalized)
}

samples_columns = {
    'sample_name': 'sample_name'    # Column in the samples file for sample names (generalized) 
}

# Load the ENA file and the samples file
ENA_file = pd.read_csv('filereport_read_run_PRJNA1158265_tsv.txt', sep='\t', encoding='utf-8')
desired_samples = pd.read_csv('samples_of_interest.txt', header=None, names=['sample_name'])

# Validate if the required columns exist in the ENA and sample files
for col in ENA_columns.values():
    if col not in ENA_file.columns:
        raise ValueError(f"Missing required column in ENA file: {col}")

for col in samples_columns.values():
    if col not in desired_samples.columns:
        raise ValueError(f"Missing required column in samples file: {col}")

# Create the base output directory if it doesn't exist
base_output_dir = "downloaded_files"
os.makedirs(base_output_dir, exist_ok=True)

# Function to submit Slurm jobs
def submit_slurm_job(sample_name, fastq_urls, library_id, sample_dir):
    for idx, url in enumerate(fastq_urls.split(';')):
        # Prepare the filename for the download
        output_file = os.path.join(sample_dir, f"{library_id}_{idx+1}.fastq")

        # If the URL indicates the file is gzipped, set the output file name with a .gz extension
        if url.endswith('.gz'):
            output_file = os.path.join(sample_dir, f"{library_id}_{idx+1}.fastq.gz")
        
        # Create the Slurm command that will run wget inside the Slurm job
        slurm_command = f"""
        sbatch --job-name=download_{sample_name}_{library_id}_{idx+1} --wrap="wget -q -O {output_file} {url}"
        """

        # Run the Slurm command to submit the job
        print(f"Submitting Slurm job: {slurm_command}")
        subprocess.run(slurm_command, shell=True)

# Loop through each sample to find matching rows and submit Slurm jobs for downloading
for _, sample_row in desired_samples.iterrows():
    sample_name = sample_row[samples_columns['sample_name']]  # Use the generalized column for sample name
    
    # Find all matching rows in the ENA file for this sample
    matching_rows = ENA_file[ENA_file[ENA_columns['ENA_file_id']].str.contains(sample_name)]
    
    if matching_rows.empty:
        print(f"No matching rows found for sample: {sample_name}")
        continue
    
    # Create a subdirectory for the sample
    sample_dir = os.path.join(base_output_dir, sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    # Process each matching row
    for _, matching_row in matching_rows.iterrows():
        fastq_urls = matching_row[ENA_columns['ftp_links']]
        library_id = matching_row[ENA_columns['ENA_file_id']]
        
        # Submit Slurm jobs for downloading files
        submit_slurm_job(sample_name, fastq_urls, library_id, sample_dir)

print("Slurm jobs have been submitted for file download.")
