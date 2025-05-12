import os
from pathlib import Path
import pandas as pd
from collections import defaultdict

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

# Process coverage files
coverage_results = []
for path in get_prioritized_files(coverage_file):
    print(f"\nProcessing coverage file: {path}")
    try:
        df = pd.read_csv(path, sep="\t")
        print(f"  Columns found: {list(df.columns)}")
        if {'Sample', 'bam_file', 'mean_coverage'}.issubset(df.columns):
            df_subset = df[['Sample', 'mean_coverage', 'bam_file']].copy()
            print(f"  Number of rows before BAM resolution: {len(df_subset)}")

            # Get base directory above multiqc_data
            base_dir = path.parents[2]

            # Resolve BAM paths
            df_subset['bam_file'] = df_subset['bam_file'].apply(
                lambda x: find_bam_path(x, base_dir)
            )

            print(f"  Number of rows after BAM resolution: {df_subset['bam_file'].notna().sum()} (non-null)")
            coverage_results.append(df_subset)
        else:
            print(f"Warning: Required columns missing in {path}")
    except Exception as e:
        print(f"Error reading {path}: {e}")

# Combine coverage
if coverage_results:
    combined_df = pd.concat(coverage_results, ignore_index=True)
    print(f"\nCombined rows before dropna: {len(combined_df)}")

    # Drop rows where bam_file could not be resolved
    combined_df = combined_df.dropna(subset=['bam_file'])
    print(f"Rows after dropna: {len(combined_df)}")

    # Drop duplicates
    before_dedup = len(combined_df)
    combined_df.drop_duplicates(subset=['Sample', 'bam_file'], inplace=True)
    print(f"Rows after deduplication: {len(combined_df)} (removed {before_dedup - len(combined_df)} duplicates)")
else:
    raise RuntimeError("No valid coverage files found.")

# Save output
combined_df.to_csv('merged_coverage_duplication.csv', index=False)
print("\nSaved to merged_coverage_duplication.csv")
