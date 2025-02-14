#!/bin/bash
##SBATCH --cpus-per-task=1 #should be number of cores you will be using
#SBATCH --time=4-00:00
##S#BATCH --job-name=

module load nextflow 

nextflow run variant_calling.nf -resume -process.echo  -c /gpfs/home/xrq24scu/fox_repo/variant_calling/variantcalling-config
    #tee process.txt # tee is not necessary if you submit this as a bash job
