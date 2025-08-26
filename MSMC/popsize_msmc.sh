#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=msmc
#SBATCH --output=msmc_%A_%a.out
#SBATCH --error=msmc_%A_%a.err


##############MODIFY THIS TO USE the output of ordered_df_filtered, after pushing to github


###submit with: sbatch --array=1-5 popsize_msmc.sh
### (replace 5 with the number of populations you have)

MSMC_OUTDIR=/gpfs/data/bergstrom/paula/fox_repo/MSMC
INPUT_DIR=/gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering/autosomes/msmc_files/multihetsepfiles

# Population number for this task (assumes SLURM_ARRAY_TASK_ID starts at 1)
POP_NUM=$SLURM_ARRAY_TASK_ID

# Get all samples belonging to this population
SAMPLES=($(awk -v pop="$POP_NUM" '$2==pop {print $1}' pop_groups.txt))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "ERROR: No samples found for population $POP_NUM"
    exit 1
fi

echo "Running MSMC for population $POP_NUM using samples: ${SAMPLES[@]}"

# Gather input files
FILES=()
for S in "${SAMPLES[@]}"; do
    FILES+=( ${INPUT_DIR}/${S}.chr*.multihetsep.txt )
done

# Check that all files exist
for F in "${FILES[@]}"; do
    if [ ! -f $F ]; then
        echo "ERROR: Input file not found: $F"
        exit 1
    fi
done

# Convert array to space-separated string for MSMC
FILES_STR="${FILES[@]}"

# Run MSMC
msmc -t 4 -o "${MSMC_OUTDIR}/pop${POP_NUM}.msmc" $FILES_STR
