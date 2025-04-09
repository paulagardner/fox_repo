#module load python/anaconda/2024.06/3.12.4
import pandas as pd
import os
import subprocess

#module load python/anaconda/2024.06/3.12.4

# Define required column mappings
ENA_columns = {
    'sample_alias': 'sample_alias',  # Column in the ENA file for sample name
    'ftp_links': 'fastq_ftp',        # Column in the ENA file for fastq URLs
    'run_accession': 'run_accession' # Column containing the base file name
}

# Load ENA and samples files
ENA_file = pd.read_csv('filereport_read_run_PRJNA951250_tsv.txt', sep='\t', encoding='utf-8')
desired_samples = pd.read_csv('samples_of_interest.txt', header=None, names=['sample_alias'])

# Convert sample_alias to strings for proper matching
desired_samples['sample_alias'] = desired_samples['sample_alias'].astype(str)
ENA_file['sample_alias'] = ENA_file['sample_alias'].astype(str)

# Validate required columns
for col in ENA_columns.values():
    if col not in ENA_file.columns:
        raise ValueError(f"Missing required column in ENA file: {col}")

# Create base output directory
base_output_dir = "downloaded_files"
os.makedirs(base_output_dir, exist_ok=True)

# Function to submit Slurm jobs
def submit_slurm_job(sample_name, fastq_urls, run_accession, sample_dir):
    fastq_files = fastq_urls.split(';')  # Split URLs for paired-end reads

    # Remove any existing .fastq or .fastq.gz extension from run_accession to form a base name
    if run_accession.endswith('.fastq.gz'):
        base_name = run_accession[:-9]  # Remove '.fastq.gz'
    elif run_accession.endswith('.fastq'):
        base_name = run_accession[:-6]  # Remove '.fastq'
    else:
        base_name = run_accession

    for idx, url in enumerate(fastq_files):
        # Choose the proper extension based on the URL
        if url.endswith('.fastq.gz'):
            ext = '.fastq.gz'
        else:
            ext = '.gz'
        
        # For paired-end files, append _1, _2, etc.; for single files, no extra numbering is added.
        if len(fastq_files) > 1:
            filename = f"{base_name}_{idx+1}{ext}"
        else:
            filename = f"{base_name}{ext}"
        
        output_file = os.path.join(sample_dir, filename)

        # Build the Slurm command.
        # Using single quotes around the wget command can help avoid shell quoting issues.
        slurm_command = (
            f"sbatch --job-name=download_{sample_name}_{idx+1} "
            f"--cpus-per-task=4 --partition=compute-64-512 "
            f"--wrap='wget -q -O {output_file} {url}'"
        )

        print(f"Submitting Slurm job: {slurm_command}")
        # Submit the job and capture output for debugging
        result = subprocess.run(slurm_command, shell=True, capture_output=True, text=True)
        print("Return code:", result.returncode)
        if result.stdout:
            print("stdout:", result.stdout)
        if result.stderr:
            print("stderr:", result.stderr)

# Process each sample
for _, sample_row in desired_samples.iterrows():
    sample_name = sample_row['sample_alias']
    
    # Find matching rows in ENA file based on sample_alias
    matching_rows = ENA_file[ENA_file['sample_alias'] == sample_name]

    if matching_rows.empty:
        print(f"No matching rows found for sample: {sample_name}")
        continue
    
    # Create sample directory (sample_alias folder)
    sample_dir = os.path.join(base_output_dir, sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    # Submit jobs for each match
    for _, matching_row in matching_rows.iterrows():
        fastq_urls = matching_row['fastq_ftp']
        run_accession = matching_row['run_accession']  # Base filename
        submit_slurm_job(sample_name, fastq_urls, run_accession, sample_dir)

print("Slurm jobs have been submitted for file download.")
