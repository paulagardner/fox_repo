import os
from pathlib import Path
import pandas as pd
from collections import defaultdict

# Base directory
base_dir = Path("/gpfs/data/bergstrom/paula/fox_repo")

# File names to look for
coverage_file = "multiqc_qualimap_bamqc_genome_results.txt"
duplication_file = "multiqc_picard_dups.txt"
qualimap_file = "genome_results.txt"

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


def get_qualimap_files(filename):
    found = defaultdict(list)

    for root, dirs, files in os.walk(base_dir, followlinks=True):
        if filename in files:
            full_path = os.path.join(root, filename)
            sample_dir = os.path.normpath(full_path.split('/qualimap')[0])
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
#for path in get_prioritized_files(coverage_file):
#    print(f"\nProcessing coverage file: {path}")

for path in get_qualimap_files(qualimap_file):
    print(f"\nProcessing coverage and dup file: {path}")
