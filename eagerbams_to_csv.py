import os
from pathlib import Path
import pandas as pd
from collections import defaultdict
import pysam  # <-- added

# Base directory
base_dir = Path("/gpfs/data/bergstrom/paula/fox_repo")

# File targets
coverage_file = "multiqc_qualimap_bamqc_genome_results.txt"
duplication_file = "multiqc_picard_dups.txt"

# Store best file paths (prioritize eager_results)
def get_prioritized_files(filename):
    found = defaultdict(list)

    for root, dirs, files in os.walk(base_dir, followlinks=True):
        if filename in files:
            full_path = os.path.join(root, filename)

            # Print out the symlink target
            print(f"Found file: {full_path}")
            if os.path.islink(full_path):
                print(f"  Symlink target: {os.readlink(full_path)}")

            sample_dir = os.path.normpath(full_path.split('/multiqc_data')[0])
            found[sample_dir].append(full_path)

    final_paths = []
    for sample_dir, paths in found.items():
        # Prefer logical paths that include 'eager_results'
        preferred = [p for p in paths if 'eager_results' in str(p)]
        chosen = preferred[0] if preferred else paths[0]
        print(f"Chosen for {sample_dir}: {chosen}")
        final_paths.append(Path(chosen))

    return final_paths


# Function to find BAM file in subdirectories below the base directory
def find_bam_path(bam_filename, base_dir):
    matches = []

    for root, dirs, files in os.walk(base_dir):
        if bam_filename in files:
            matches.append(os.path.join(root, bam_filename))

    if not matches:
        print(f"Warning: BAM file {bam_filename} not found under {base_dir}")
        return None

    # Prefer matches NOT in '/work/' directories
    matches.sort(key=lambda x: '/work/' in x)
    print(f"Resolved {bam_filename} to {matches[0]}")
    return matches[0]


# Extract sample name from BAM read group
def get_sample_name_from_bam(bam_path):
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        read_groups = bam.header.get('RG', [])
        if not read_groups:
            print(f"Warning: No read group (RG) found in BAM header: {bam_path}")
            return None
        # Get first RG SM tag
        sample_name = read_groups[0].get('SM')
        if sample_name:
            return sample_name
        else:
            print(f"Warning: No SM tag in read group for BAM: {bam_path}")
            return None
    except Exception as e:
        print(f"Error reading BAM file {bam_path}: {e}")
        return None


# Process coverage files
coverage_results = []
for path in get_prioritized_files(coverage_file):
    print(f"\nProcessing coverage file: {path}")
    try:
        df = pd.read_csv(path, sep="\t")
        print(f"  Columns found: {list(df.columns)}")
        if {'Sample', 'bam_file', 'mean_coverage'}.issubset(df.columns):
            df_subset = df[['mean_coverage', 'bam_file']].copy()
            print(f"  Number of rows before BAM resolution: {len(df_subset)}")

            # Get base directory above multiqc_data
            base_dir = path.parents[2]

            # Resolve BAM paths
            df_subset['bam_file'] = df_subset['bam_file'].apply(
                lambda x: find_bam_path(x, base_dir)
            )

            # Drop rows without resolved BAMs
            df_subset = df_subset.dropna(subset=['bam_file'])

            # Extract sample name from BAM
            df_subset['Sample'] = df_subset['bam_file'].apply(get_sample_name_from_bam)
            df_subset = df_subset.dropna(subset=['Sample'])

            print(f"  Number of rows after BAM + sample resolution: {len(df_subset)}")
            coverage_results.append(df_subset)
        else:
            print(f"Warning: Required columns missing in {path}")
    except Exception as e:
        print(f"Error reading {path}: {e}")

# Combine coverage
if coverage_results:
    combined_df = pd.concat(coverage_results, ignore_index=True)
    print(f"\nCombined rows before deduplication: {len(combined_df)}")

    # Reorder columns as Sample, bam_file, mean_coverage
    combined_df = combined_df[['Sample', 'bam_file', 'mean_coverage']]


    # Drop duplicates
    before_dedup = len(combined_df)
    combined_df.drop_duplicates(subset=['Sample', 'bam_file'], inplace=True)
    print(f"Rows after deduplication: {len(combined_df)} (removed {before_dedup - len(combined_df)} duplicates)")
else:
    raise RuntimeError("No valid coverage files found.")

# Save output
combined_df.to_csv('merged_coverage_duplication.csv', index=False)
print("\nSaved to merged_coverage_duplication.csv")