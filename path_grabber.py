#!/usr/bin/env python/anaconda/2024.06/3.12.4
# module load python/anaconda/2024.06/3.12.4
import csv
from pathlib import Path

def extract_sample_id(bam_filename):
    """
    Extracts the sample id from a bam filename.
    Removes everything from '_rmdup' onward, then splits the remainder on underscores.
    If there are at least two tokens, returns the first two joined by an underscore.
    """
    base = bam_filename.split('_rmdup')[0]
    parts = base.split('_')
    if len(parts) >= 2:
        return "_".join(parts[:2])
    else:
        return base

def extract_sample_from_multiqc(sample_str):
    """
    Extracts a sample id from the multiqc 'Sample' column by taking the first two underscore-separated tokens.
    """
    parts = sample_str.split('_')
    if len(parts) >= 2:
        return "_".join(parts[:2])
    else:
        return sample_str

def get_group_dir(bam_path: Path) -> str:
    """
    Returns the directory name immediately preceding 'eager_results' in the BAM file's path.
    If 'eager_results' is not found, returns the immediate parent directory's name.
    """
    parts = bam_path.parts
    if "eager_results" in parts:
        idx = parts.index("eager_results")
        if idx > 0:
            return parts[idx-1]
    return bam_path.parent.name

# 1. Set the base directory (common root) containing both BAMs and multiqc files.
base_dir = Path("/gpfs/home/xrq24scu/fox_repo")

# 2. Gather coverage data from all multiqc output files.
coverage_by_sample = {}
# Search recursively for any multiqc file matching the pattern.
gr_files = list(base_dir.rglob("**/eager_results/multiqc/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt"))
print("Found multiqc files:", len(gr_files))

for gr_file in gr_files:
    print("Processing multiqc file:", gr_file)
    with open(gr_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        try:
            sample_index = header.index("Sample")
            coverage_index = header.index("mean_coverage")
        except ValueError as e:
            print("Header error in", gr_file, ":", e)
            continue
        for row in reader:
            if not row or len(row) <= max(sample_index, coverage_index):
                continue
            raw_sample = row[sample_index].strip()
            sample_id = extract_sample_from_multiqc(raw_sample)
            coverage = row[coverage_index].strip()
            coverage_by_sample[sample_id] = coverage
            print("MultiQC: sample", sample_id, "coverage", coverage)

print("Final coverage dictionary:", coverage_by_sample)

# 3. Gather BAM files from both merged_bams and deduplication folders.
#    If a sample appears in both locations, only the merged_bams file is used.
bam_dict = {}

# Priority: merged_bams
merged_bam_files = list(base_dir.rglob("eager_results/merged_bams/**/*_rmdup.bam"))
print("Found merged_bams files:", len(merged_bam_files))
for bam_file in merged_bam_files:
    sample_id = extract_sample_id(bam_file.name)
    if sample_id:
        bam_dict[sample_id] = bam_file

# Then add deduplication files only if the sample isn’t already present.
dedup_bam_files = list(base_dir.rglob("**/deduplication/**/*_rmdup.bam"))
print("Found deduplication files:", len(dedup_bam_files))
for bam_file in dedup_bam_files:
    sample_id = extract_sample_id(bam_file.name)
    if sample_id and sample_id not in bam_dict:
        bam_dict[sample_id] = bam_file

# 4. Write output CSV with: sample, path, mean_coverage.
#    Sorting output by the directory immediately above eager_results (group) then by filename.
output_csv = "all_bam_path_info.csv"
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample", "path", "mean_coverage"])
    
    # Create a sorted list of items.
    sorted_items = sorted(
        bam_dict.items(),
        key=lambda item: (
            get_group_dir(item[1]).lower(),
            item[1].name.lower()
        )
    )
    
    # Debug: print group for each file
    for sample_id, bam_file in sorted_items:
        group = get_group_dir(bam_file)
        print(f"Sample {sample_id} in group '{group}' from file {bam_file.name}")
    
    for sample_id, bam_file in sorted_items:
        # Try the extracted sample_id from the two-token rule.
        mean_coverage = coverage_by_sample.get(sample_id)
        if not mean_coverage or mean_coverage == "":
            # Fallback: try using only the first token
            alt_sample = bam_file.name.split('_')[0]
            mean_coverage = coverage_by_sample.get(alt_sample, "NA")
        writer.writerow([sample_id, str(bam_file), mean_coverage])
        print("Output:", sample_id, str(bam_file), mean_coverage)
