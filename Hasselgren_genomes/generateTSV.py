import os
import pandas as pd
import subprocess

#module load python/anaconda/2024.06/3.12.4

# Load the ENA file
ENA_file = pd.read_csv('filereport_read_run_PRJEB76449_tsv.txt', sep='\t', encoding='utf-8')

# Define columns
ENA_columns = {
    'run_accession': 'run_accession',
    'scientific_name': 'scientific_name',
    'instrument_model': 'instrument_model'
}

download_dir = "downloaded_files"

def get_lane_number(file_path):
    """Extract Lane number using zcat only if the file exists."""
    if not os.path.exists(file_path):
        return "NA"
    try:
        output = subprocess.check_output(f"zcat {file_path} | head -1", shell=True, text=True)
        return output.split(':')[3]  # Extract the single digit between two colons
    except Exception as e:
        print(f"Error extracting lane number from {file_path}: {e}")
        return "NA"

# Find all downloaded .fastq.gz files
file_mapping = {}
for root, _, files in os.walk(download_dir):
    for file in files:
        if file.endswith(".fastq.gz"):
            full_path = os.path.join(root, file)
            run_accession = file.split("_")[0]  # Extract run_accession from filename
            sample_name = os.path.basename(os.path.dirname(full_path))  # Extract immediate parent directory as Sample_Name

            if run_accession not in file_mapping:
                file_mapping[run_accession] = {'files': [], 'sample_name': sample_name}
            
            file_mapping[run_accession]['files'].append(full_path)

# Prepare output data
output_data = []
for run_accession, data in file_mapping.items():
    file_paths = data['files']
    sample_name = data['sample_name']  # Assign Sample_Name based on parent directory

    # Find matching row in ENA file
    matching_row = ENA_file[ENA_file[ENA_columns['run_accession']].astype(str) == str(run_accession)]
    organism = matching_row[ENA_columns['scientific_name']].values[0] if not matching_row.empty else "Unknown"
    instrument_model = matching_row[ENA_columns['instrument_model']].values[0] if not matching_row.empty else "Unknown"

    # Assign Colour_Chemistry based on instrument model
    if "NovaSeq 6000" in instrument_model:
        colour_chemistry = 2
    elif "MiSeq" in instrument_model:
        colour_chemistry = 4
    else:
        colour_chemistry = "NA"

    r1_file = next((f for f in file_paths if "_1" in f), "NA")
    r2_file = next((f for f in file_paths if "_2" in f), "NA")
    lane = get_lane_number(r1_file if r1_file != "NA" else r2_file)

    output_data.append({
        'Sample_Name': sample_name,  # Now using extracted parent directory name
        'Library_ID': sample_name,  # Library ID remains the run_accession
        'Lane': lane,
        'Colour_Chemistry': colour_chemistry,
        'SeqType': 'PE',
        'Organism': organism,
        'Strandedness': 'double',
        'UDG_Treatment': 'none',
        'R1': r1_file,
        'R2': r2_file,
        'BAM': 'NA'
    })

# Create DataFrame and save as .tsv
output_df = pd.DataFrame(output_data)
output_df.to_csv('inputFile.tsv', sep='\t', index=False)
print("inputFile.tsv has been generated based on downloaded files.")
