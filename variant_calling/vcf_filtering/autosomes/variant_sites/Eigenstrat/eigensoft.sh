#!/bin/bash

module load EIGENSOFT/8.0.0
module load gsl #need to load these to get smartpca to work
module load openblas 
module load plink/2.00

plink2 --vcf /gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering/autosomes/variant_sites/Vvulpes.c2.samplesites01.vcf.gz \
        --out Vvulpes_plink_autosomal --keep-allele-order --const-fid 0

smartpca -p par.smartpca

/gpfs/data/bergstrom/sw/bin/fix-smartpca-evec-table.sh Vvulpes.pca.eval > Vvulpestidy.evec


