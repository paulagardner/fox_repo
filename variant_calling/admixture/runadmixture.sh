#!/bin/bash

# Input and output PLINK dataset prefixes (no extensions)
INPUT_PREFIX="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpes_plink_autosomal"
OUTPUT_DIR="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture"
OUTPUT_PREFIX="${OUTPUT_DIR}/low_missingness_rm_outlier"

# Threshold for missingness (e.g., remove samples with F_MISS > 0.9)
MISSINGNESS_THRESHOLD=0.9

################### Step 1: Generate Missingness Report
echo "Running PLINK to calculate missingness..."

module load plink/2.00

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Generate the missingness report
plink2 --bfile "$INPUT_PREFIX" --missing --out "${OUTPUT_PREFIX}"

# Check if the .smiss file was generated
SMISS_FILE="${OUTPUT_PREFIX}.smiss"
if [[ ! -f "$SMISS_FILE" ]]; then
  echo "ERROR: Missingness file not found at $SMISS_FILE"
  exit 1
fi

echo "Missingness file generated at: $SMISS_FILE"

################### Step 2: Extract samples to remove based on missingness threshold
TMP_REMOVE_FILE=$(mktemp)
echo "Creating temporary remove file: $TMP_REMOVE_FILE"

while IFS=$'\t' read -r FID IID MISSING_CT OBS_CT F_MISS; do
    if [[ "$FID" == "#FID" || "$FID" == "FID" ]]; then
        continue
    fi
    if (( $(echo "$F_MISS > $MISSINGNESS_THRESHOLD" | bc -l) )); then
        echo -e "$FID\t$IID" >> "$TMP_REMOVE_FILE"
        echo "Added $IID ($FID) to remove list (F_MISS=$F_MISS)"
    fi
done < "$SMISS_FILE"

# Manually add a specific sample
echo -e "0\tYPI1082" >> "$TMP_REMOVE_FILE"
echo "Manually added YPI1082 to remove list"

if [[ ! -s "$TMP_REMOVE_FILE" ]]; then
    echo "No samples exceed the missingness threshold. Nothing to remove."
else
    echo "Samples to be removed:"
    cat "$TMP_REMOVE_FILE"
fi

echo "Contents of remove file:"
cat "$TMP_REMOVE_FILE"

################### Step 3: Remove samples
echo "Running PLINK to remove specified samples..."

plink2 --bfile "$INPUT_PREFIX" \
       --remove "$TMP_REMOVE_FILE" \
       --make-bed \
       --out "$OUTPUT_PREFIX" \
       --allow-no-sex

echo "PLINK finished. Filtered dataset saved as: $OUTPUT_PREFIX"
echo "Remaining samples after removal:"
plink2 --bfile "$OUTPUT_PREFIX" --freq

################### Step 3.5: Output sample lists
REMOVED_SAMPLES_FILE="${OUTPUT_PREFIX}_removed_samples.txt"
KEPT_SAMPLES_FILE="${OUTPUT_PREFIX}_final_samples.txt"

cut -f2 "$TMP_REMOVE_FILE" | sort > "$REMOVED_SAMPLES_FILE"
echo "Removed samples list written to: $REMOVED_SAMPLES_FILE"

cut -f2 "${OUTPUT_PREFIX}.fam" | sort > "$KEPT_SAMPLES_FILE"
echo "Final retained samples list written to: $KEPT_SAMPLES_FILE"

# Clean up
echo "Cleaning up temporary remove file..."
rm "$TMP_REMOVE_FILE"
echo "Temporary file cleaned up."



##############filter genotypes that are completely missing due to having only 
##been present, presumably, in the YPI1082 dataset

plink2 --bfile "$OUTPUT_PREFIX" --geno 0.999 --make-bed --out "${OUTPUT_PREFIX}_snpmanualfilter"

################### Step 4: Run ADMIXTURE
echo "Submitting ADMIXTURE jobs for K values from 3 to 10..."
echo "Using input: ${OUTPUT_PREFIX}.bed"

for K in {3..11}; do
  echo "Submitting ADMIXTURE job for K=$K"
  sbatch --job-name=admix_K$K \
         --output="${OUTPUT_DIR}/log${K}.out" \
         --time=2-1:00:00 \
         --cpus-per-task=4 \
         --mem=45G \
         --wrap="admixture --cv ${OUTPUT_PREFIX}_snpmanualfilter.bed $K"
done

echo "All ADMIXTURE jobs submitted."
