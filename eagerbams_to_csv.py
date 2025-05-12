import os
from pathlib import Path
import pandas as pd
from collections import defaultdict

# Base directory
base_dir = Path("/gpfs/data/bergstrom/paula/fox_repo")

# File names to look for
coverage_file = "multiqc_qualimap_bamqc_genome_results.txt"
duplication_file = "multiqc_picard_dups.txt"

# Containers for results
coverage_results = []
duplication_results = []

# Get prioritized file paths
def get_prioritized_files(filename):
    found = defaultdict(list)

    for root, dirs, files in os.walk(base_dir, followlinks=True):
        if filename in files:
            full_path = os.path.join(root, filename)
            sample_dir = os.path.normpath(full_path.split('/multiqc_data')[0])
            found[sample_dir].append(full_path)

    final_paths = []
    for sample_dir, paths in found.items():
        preferred = [p for p in paths if 'eager_results' in str(p)]
        chosen = preferred[0] if preferred else paths[0]
        print(f"Chosen for {sample_dir}: {chosen}")
        final_paths.append(Path(chosen))

    return final_paths

# Resolve BAM paths
def find_bam_path(bam_filename, base_dir):
    matches = []
    for root, dirs, files in os.walk(base_dir):
        if bam_filename in files:
            matches.append(os.path.join(root, bam_filename))

    if not matches:
        print(f"Warning: BAM file {bam_filename} not found under {base_dir}")
        return None

    matches.sort(key=lambda x: '/work/' in x)
    print(f"Resolved {bam_filename} to {matches[0]}")
    return matches[0]

# Process each coverage file
for path in get_prioritized_files(coverage_file):
    print(f"\nProcessing coverage file: {path}")
    try:
        df = pd.read_csv(path, sep="\t")
        if {'Sample', 'bam_file', 'mean_coverage'}.issubset(df.columns):
            df_subset = df[['Sample', 'mean_coverage', 'bam_file']].copy()
            base_path = path.parents[2]  # Get directory above multiqc_data

            df_subset['bam_file'] = df_subset['bam_file'].apply(
                lambda x: find_bam_path(x, base_path)
            )
            df_subset = df_subset.dropna(subset=['bam_file'])
            coverage_results.append(df_subset)

            # Now check for a duplication file in the same dir
            dup_path = path.parent / duplication_file
            if dup_path.exists():
                print(f"Processing duplication file: {dup_path}")
                try:
                    dup_df = pd.read_csv(dup_path, sep="\t")
                    dup_df.columns = [c.upper() for c in dup_df.columns]
                    if {'SAMPLE', 'PERCENT_DUPLICATION'}.issubset(dup_df.columns):
                        dup_df = dup_df[['SAMPLE', 'PERCENT_DUPLICATION']].copy()
                        dup_df = dup_df.rename(columns={'SAMPLE': 'Sample', 'PERCENT_DUPLICATION': 'percent_duplication'})
                        duplication_results.append(dup_df)
                except Exception as e:
                    print(f"Error reading {dup_path}: {e}")
            else:
                print(f"No duplication file found at {dup_path}")

        else:
            print(f"Warning: Required columns missing in {path}")
    except Exception as e:
        print(f"Error reading {path}: {e}")

# Combine coverage data
if not coverage_results:
    raise RuntimeError("No valid coverage files found.")

combined_df = pd.concat(coverage_results, ignore_index=True)
print(f"\nCombined rows before deduplication: {len(combined_df)}")
combined_df.drop_duplicates(subset=['Sample', 'bam_file'], inplace=True)
print(f"Rows after deduplication: {len(combined_df)}")

# Merge duplication if available
if duplication_results:
    duplication_df = pd.concat(duplication_results, ignore_index=True)
    duplication_df['Sample'] = duplication_df['Sample'].astype(str)
    combined_df['Sample'] = combined_df['Sample'].astype(str)

    combined_df = pd.merge(combined_df, duplication_df, on='Sample', how='left')
    print(f"Combined duplication rows: {len(duplication_df)}")
else:
    print("No valid duplication data found. Proceeding without it.")

# Reorder and save
combined_df = combined_df[['Sample', 'bam_file', 'mean_coverage', 'percent_duplication']]
combined_df.to_csv('merged_coverage_duplication.csv', index=False)
print("\nSaved to merged_coverage_duplication.csv")
