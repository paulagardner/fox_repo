#!/bin/bash

module load EIGENSOFT/8.0.0
module load gsl #need to load these to get smartpca to work
module load openblas 
module load plink/2.00


########## I should filetr for missingness much earlier- at this step- so it doesn't have to be done later and then samples
### with too much missing data removed from the Rscript for PCA via Eigensoft results, going backwards ina  confusing way

plink2 --vcf /gpfs/data/bergstrom/paula/fox_repo/variant_calling/vcf_filtering/autosomes/variant_sites/Vvulpes.c2.samplesites01.vcf.gz \
        --make-bed --out Vvulpes_plink_autosomal --keep-allele-order --const-fid 0

smartpca -p /gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/par.smartpca

/gpfs/data/bergstrom/sw/bin/fix-smartpca-evec-table.sh Vvulpes.pca.evec > Vvulpestidy.evec

#########repeat for no-outliers
smartpca -p /gpfs/data/bergstrom/paula/fox_repo/variant_calling/Eigenstrat/par_removeoutliers.smartpca

/gpfs/data/bergstrom/sw/bin/fix-smartpca-evec-table.sh Vvulpes_rmoutliers.pca.evec > Vvulpes_rmoutliers_tidy.evec


