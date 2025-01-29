import os
import pandas as pd
import subprocess

##first, start with acivating conda via module load python/anaconda/2024.06/3.12.4

# Load the ENA file
ENA_file = pd.read_csv('filereport_read_run_PRJNA1158265_tsv.txt', sep='\t', encoding='utf-8')

# Define columns
ENA_columns = {
    'ENA_file_id': 'library_name',
    'scientific_name': 'scientific_name'  # Organism name
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
            library_id = file.split("_")[0]  # Extract library_id from filename
            sample_id = os.path.basename(os.path.dirname(full_path))  # Directory above the file
            if library_id not in file_mapping:
                file_mapping[library_id] = {'sample_id': sample_id, 'files': []}
            file_mapping[library_id]['files'].append(full_path)

# Prepare output data
output_data = []
for library_id, data in file_mapping.items():
    file_paths = data['files']
    sample_id = data['sample_id']
    matching_row = ENA_file[ENA_file[ENA_columns['ENA_file_id']] == library_id]
    organism = matching_row[ENA_columns['scientific_name']].values[0] if not matching_row.empty else "Unknown"
    
    r1_file = next((f for f in file_paths if "_1" in f), "NA")
    r2_file = next((f for f in file_paths if "_2" in f), "NA")
    lane = get_lane_number(r1_file if r1_file != "NA" else r2_file)
    
    output_data.append({
        'Sample_Name': sample_id,
        'Library_ID': library_id,
        'Lane': lane,
        'Colour_Chemistry': 2,
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
