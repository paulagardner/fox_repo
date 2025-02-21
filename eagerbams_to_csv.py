#!/usr/bin/env python/anaconda/2024.06/3.12.4
# module load python/anaconda/2024.06/3.12.4
import csv
from pathlib import Path

def extract_sample_id(bam_filename):
    """
    Extracts the sample id from a BAM filename.
    Removes everything from '_rmdup' onward, then splits on underscores.
    If the second token equals "udgnone" (case-insensitive), returns only the first token.
    Otherwise, if there are at least two tokens, returns the first two joined by an underscore.
    """
    base = bam_filename.split('_rmdup')[0]
    parts = base.split('_')
    if len(parts) >= 2:
        if parts[1].lower() == "udgnone":
            return parts[0]
        else:
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

def normalize_sample(s):
    """
    Normalizes a sample identifier by removing non-alphanumeric characters.
    """
    return "".join(ch for ch in s if ch.isalnum())

def get_group_from_path(path: Path) -> str:
    """
    Returns a group tag based on the path.
    For example, if 'our_genomes' is in the path, returns "our_genomes";
    if 'NA_genomes' is in the path, returns "NA_genomes"; otherwise, returns "unknown".
    """
    s = path.as_posix()
    if "our_genomes" in s:
        return "our_genomes"
    elif "NA_genomes" in s:
        return "NA_genomes"
    else:
        return "unknown"

# 1. Set the base directory containing BAMs, multiqc files, and inputFile.tsv files.
base_dir = Path("/gpfs/home/xrq24scu/fox_repo")

# 2. Gather coverage data from all multiqc output files.
coverage_by_sample = {}
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

# 3. Build library mapping from inputFile.tsv files, tagging each mapping with its group.
library_mapping = []
input_tsv_files = list(base_dir.rglob("inputFile.tsv"))
print("Found inputFile.tsv files:", len(input_tsv_files))
for tsv_file in input_tsv_files:
    mapping_group = get_group_from_path(tsv_file)
    print("Processing mapping file:", tsv_file, "group:", mapping_group)
    with open(tsv_file, "r") as tsv_in:
        reader = csv.DictReader(tsv_in, delimiter="\t")
        for row in reader:
            lib_id = row.get("Library_ID")
            sample_name = row.get("Sample_Name")
            if lib_id and sample_name:
                # Store a tuple (Library_ID, Sample_Name, group)
                library_mapping.append((lib_id, sample_name, mapping_group))
                print(f"Mapping: Library_ID '{lib_id}' -> Sample_Name '{sample_name}' (group: {mapping_group})")
            else:
                print("Row missing Library_ID or Sample_Name:", row)
print(f"Total mapping entries collected: {len(library_mapping)}")

# 4. Gather BAM files from merged_bams and deduplication folders.
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

# 5. Write output CSV with: sample, path, mean_coverage.
#    Sorted by the directory immediately above eager_results (group) then by filename.
output_csv = "all_bam_path_info.csv"
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample", "path", "mean_coverage"])
    
    sorted_items = sorted(
        bam_dict.items(),
        key=lambda item: (
            get_group_dir(item[1]).lower(),
            item[1].name.lower()
        )
    )
    
    for sample_id, bam_file in sorted_items:
        # Lookup coverage using the raw sample ID.
        mean_coverage = coverage_by_sample.get(sample_id)
        if not mean_coverage or mean_coverage == "":
            alt_sample = bam_file.name.split('_')[0]
            mean_coverage = coverage_by_sample.get(alt_sample, "NA")
        
        # Determine the BAM file's group.
        bam_group = get_group_from_path(bam_file)
        
        # Apply library mapping using substring match, but only consider mappings from the same group.
        final_sample = sample_id
        norm_sample = normalize_sample(sample_id)
        for lib_id, mapped_name, map_group in library_mapping:
            if map_group == bam_group and norm_sample in normalize_sample(lib_id):
                final_sample = mapped_name
                print(f"Mapping: replacing raw sample '{sample_id}' with '{mapped_name}' (group: {map_group})")
                break
        
        writer.writerow([final_sample, str(bam_file), mean_coverage])
        print("Output:", final_sample, str(bam_file), mean_coverage)
