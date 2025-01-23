#!/bin/bash
##SBATCH --cpus-per-task=1 #should be number of cores you will be using
##S#BATCH --job-name=myworkflow


module load nextflow 

nextflow run bamfile_info.nf -resume -process.echo #\
    #tee process.txt # tee is not necessary if you submit this as a bash job
