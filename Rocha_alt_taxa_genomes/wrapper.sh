#!/bin/bash
#SBATCH --cpus-per-task=1 #should be number of cores you will be using
#SBATCH --time=7-00:00
#SBATCH --job-name=eager


# if you run nextflow from some environmental module or conda env, load it here.
module load nextflow/22.10.6

module load singularity 

nextflow run nf-core/eager \
    -with-singularity \
    -c /gpfs/home/xrq24scu/fox_repo/our_genomes/eager_mapping/ada-config \
    -resume \
    --publish_dir_mode 'symlink' \
    --input inputFile.tsv \
    --outdir /gpfs/data/bergstrom/paula/fox_repo/Rocha_alt_taxa_genomes/eager_results \
    --fasta /gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa \
    --fasta_index /gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa.fai \
    --seq_dict /gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.dict \
    --bwa_index /gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/ \
    --mapper bwamem \
    --skip_collapse \
    --skip_damage_calculation \
    --skip_preseq \
