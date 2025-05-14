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


import re

def parse_qualimap_file(filepath):
    """Extract sample name, bam file name, mean coverage, and duplication rate from a genome_results.txt"""
    with open(filepath, 'r') as f:
        text = f.read()

    # Get bam file name
    bam_match = re.search(r"bam file\s+=\s+(.+)", text)
    bam_filename = bam_match.group(1).strip() if bam_match else None

    # Get mean coverage
    coverage_match = re.search(r"mean coverageData\s*=\s*([\d\.]+)X", text)
    mean_coverage = float(coverage_match.group(1)) if coverage_match else None

    # Get duplication rate
    dup_match = re.search(r"duplication rate\s*=\s*([\d\.]+)%", text)
    percent_duplicates = float(dup_match.group(1)) if dup_match else None

    return bam_filename, mean_coverage, percent_duplicates


def find_bam_path_nearby(sample_dir, bam_filename):
    """Search for a BAM file in the sample_dir and its subdirectories"""
    for root, dirs, files in os.walk(sample_dir):
        if bam_filename in files:
            return os.path.join(root, bam_filename)
    print(f"Warning: BAM file {bam_filename} not found in {sample_dir}")
    return None


# Store parsed data
records = []

for qpath in get_qualimap_files("genome_results.txt"):
    print(f"Processing: {qpath}")

    # Get sample directory (i.e., the one above eager_results)
    sample_dir = Path(qpath).parts
    if "eager_results" in sample_dir:
        sample_root = Path(*sample_dir[:sample_dir.index("eager_results")])
    else:
        continue  # Skip unexpected structure

    sample_name = Path(qpath).parent.name.split("_")[0]  # e.g. "AK_S12-1159"
    bam_file, coverage, duplication = parse_qualimap_file(qpath)

    if not bam_file:
        print(f"Skipping {qpath} — no BAM info found")
        continue

    bam_path = find_bam_path_nearby(sample_root, Path(bam_file).name)

    records.append({
        "sample": sample_name,
        "bam_path": bam_path,
        "mean_coverage": coverage,
        "percent_duplicates": duplication
    })

# Convert to DataFrame and write CSV
df = pd.DataFrame(records)
df.to_csv("qualimap_summary.csv", index=False)
print("Written to qualimap_summary.csv")
