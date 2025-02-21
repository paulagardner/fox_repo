import csv
import os
import glob

# --- Parameters ---
sample_csv = "/gpfs/data/bergstrom/paula/fox_repo/all_bam_path_info.csv"  # Your sample CSV file
base_dir = "/gpfs/home/xrq24scu/fox_repo"    # Base directory to search for inputFile.tsv files
output_csv = "/gpfs/data/bergstrom/paula/fox_repo/updated_samples.csv"  # Output file for the updated CSV

print("Starting directory scan in base directory:", base_dir)

# --- Step 2: Search for inputFile.tsv files using glob ---
pattern = os.path.join(base_dir, "**", "inputFile.tsv")
print("Using glob with pattern:", pattern)
input_tsv_files = glob.glob(pattern, recursive=True)

if not input_tsv_files:
    print("No inputFile.tsv files found using glob!")
else:
    print(f"Total inputFile.tsv files found using glob: {len(input_tsv_files)}")
    for f in input_tsv_files:
        print("Found:", f)

# --- Step 3: Build a mapping from Library_ID to Sample_Name ---
library_mapping = []
for tsv_file in input_tsv_files:
    print("Processing file:", tsv_file)
    with open(tsv_file, "r") as tsv_in:
        reader = csv.DictReader(tsv_in, delimiter="\t")
        for row in reader:
            lib_id = row.get("Library_ID")
            sample_name = row.get("Sample_Name")
            if lib_id and sample_name:
                library_mapping.append((lib_id, sample_name))
                print(f"Mapping: Library_ID '{lib_id}' -> Sample_Name '{sample_name}'")
            else:
                print("Row missing Library_ID or Sample_Name:", row)

print(f"Total mapping entries collected: {len(library_mapping)}")

# --- Helper function for normalization ---
def normalize_sample(s):
    """
    Returns a normalized version of a sample identifier by removing non-alphanumeric characters.
    """
    return "".join(ch for ch in s if ch.isalnum())

# --- Step 1: Read the sample CSV file ---
print("Reading sample CSV file:", sample_csv)
with open(sample_csv, newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    sample_rows = list(reader)
print(f"Total sample rows: {len(sample_rows)}")

# --- Step 4: Update the sample names in the CSV using the mapping ---
for row in sample_rows:
    current_sample = row["sample"]
    norm_sample = normalize_sample(current_sample)
    candidate = None
    for lib_id, sample_name in library_mapping:
        # Check if the Library_ID is a substring of the current CSV sample value.
       if current_sample in lib_id:
            # If the candidate's Sample_Name (normalized) exactly matches the CSV sample, use it.
            if normalize_sample(sample_name) == norm_sample:
                candidate = sample_name
                print(f"Exact match for sample '{current_sample}' found using Library_ID '{lib_id}' with Sample_Name '{sample_name}'")
                break
            # Otherwise, if no candidate is chosen yet, store this one.
            elif candidate is None:
                candidate = sample_name
                print(f"Candidate match for sample '{current_sample}' found using Library_ID '{lib_id}' with Sample_Name '{sample_name}'")
    if candidate:
        print(f"Replacing sample '{current_sample}' with '{candidate}'")
        row["sample"] = candidate
    else:
        print(f"No mapping match for sample '{current_sample}'")

# --- Step 5: Write the updated rows back to a new CSV file ---
with open(output_csv, "w", newline="") as csv_out:
    fieldnames = ["sample", "path", "mean_coverage"]
    writer = csv.DictWriter(csv_out, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(sample_rows)