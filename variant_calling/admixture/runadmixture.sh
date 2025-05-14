#!/bin/bash

# Input and output PLINK dataset prefixes (no extensions)
INPUT_PREFIX="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/Vvulpes_plink_autosomal"
OUTPUT_PREFIX="/gpfs/data/bergstrom/paula/fox_repo/variant_calling/admixture/low_missingness_rm_outlier"

# Threshold for missingness (e.g., remove samples with F_MISS > 0.91
MISSINGNESS_THRESHOLD=0.9
################### Step 1: Generate Missingness Report
echo "Running PLINK to calculate missingness..."

module load plink/2.00
# Generate the missingness report (smiss file)
plink2 --bfile "$INPUT_PREFIX" --missing --out "$INPUT_PREFIX_missingness"

# Check if the smiss file was generated
if [[ ! -f "${INPUT_PREFIX_missingness}.smiss" ]]; then
  echo "Missingness file not found! Please check PLINK command and dataset."
  exit 1
fi

################### Step 2: Extract samples to remove based on missingness threshold
# Create a temporary file to store the list of samples to remove
TMP_REMOVE_FILE=$(mktemp)
echo "Creating temporary remove file: $TMP_REMOVE_FILE"

# Read the smiss file (assuming it is tab-delimited and no header)
while IFS=$'\t' read -r FID IID MISSING_CT OBS_CT F_MISS; do
    # Skip header line
    if [[ "$FID" == "#FID" || "$FID" == "FID" ]]; then
        continue
    fi
    # If missingness exceeds the threshold, add to the remove list
    # Using bc for floating-point comparison
    if (( $(echo "$F_MISS > $MISSINGNESS_THRESHOLD" | bc -l) )); then
        echo -e "$FID\t$IID" >> "$TMP_REMOVE_FILE"  # Correcting format to include FID
        echo "Added $IID ($FID) to remove list (F_MISS=$F_MISS)"
    fi
done < "${INPUT_PREFIX_missingness}.smiss"  # Make sure to specify the correct path to your smiss file


# Manually add a specific sample to the remove list
echo -e "YPI1082\tYPI1082" >> "$TMP_REMOVE_FILE"
echo "Manually added YPI1082 to remove list"

# Check if there are any samples to remove
if [[ ! -s "$TMP_REMOVE_FILE" ]]; then
    echo "No samples exceed the missingness threshold. Nothing to remove."
else
    echo "Samples to be removed:"
    cat "$TMP_REMOVE_FILE"
fi


################### Step 3: Run PLINK to remove samples
echo "Running PLINK to remove specified samples from the dataset..."

# Run PLINK to remove samples and write list of retained individuals
plink2 --bfile "$INPUT_PREFIX" \
       --remove "$TMP_REMOVE_FILE" \
       --make-bed \
       --out "$OUTPUT_PREFIX" \
       --allow-no-sex

# Print confirmation and the number of remaining samples after removal
echo "PLINK finished. Filtered dataset saved as: $OUTPUT_PREFIX"
echo "Remaining samples after removal:"
plink2 --bfile "$OUTPUT_PREFIX" --freq

# Clean up temporary file
echo "Cleaning up temporary remove file..."
rm "$TMP_REMOVE_FILE"
echo "Temporary file cleaned up."

################### Step 4: Run ADMIXTURE with SLURM array via sbatch --wrap

echo "Submitting ADMIXTURE jobs for K values from 3 to 11..."

for K in {3..10}; do
  echo "Submitting ADMIXTURE job for K=$K"
  sbatch --job-name=admix_K$K \
         --output=log${K}.out \
         --time=1-12:00:00 \
         --cpus-per-task=4 \
         --mem=16G \
         --wrap="admixture --cv ${OUTPUT_PREFIX}.bed $K"
done

echo "All ADMIXTURE jobs submitted."
