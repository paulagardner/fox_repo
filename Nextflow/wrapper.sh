#!/bin/bash
#SBATCH --cpus-per-task=1 #should be number of cores you will be using
#SBATCH --job-name=myworkflow


# if you run nextflow from some environmental module or conda env, load it here.
module load nextflow 
# or source /path/to/conda/bin/activate /path/to/your/env

# you might have something like $TMPDIR or $SCRATCH on your worker node - if so use it as work-dir
# if it points to /tmp you might want to create a uniquely named subdirectory
# TMPDIR=$(mktemp -d)

#nextflow run main.nf -w $TMPDIR/work -process.echo
#rm -rf $TMPDIR/work
#nextflow run test_main_processes.nf -resume -process.echo 
nextflow run read_mapping.nf --resume --process.echo
