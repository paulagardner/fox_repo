##### ORDER TO USE THESE SCRIPTS:
-----
alignment:
Ingredients needed: indexed reference genome, read groups config file containing information on the fastas.

- First: align raw fasta sequences to reference with [combined_alignment_script_nextflow_prep.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/combined_alignment_script_nextflow_prep.sh)
- sort files with [sort_bamfiles.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/sort_bamfiles.sh)
- merge files back together (as they were split up for processing) with [merge.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/merge.sh)
- mark duplicates with picard: [mark_duplicates.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/mark_duplicates.sh)
- index the merged bams with [mark_duplicates.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/index.sh)

----
post-processing:
- check that the bamfiles worked well: [validate_bamfiles.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/validate_bamfiles.sh)
- calculate coverage with picard wgs: [calculate_coverage.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/calculate_coverage.sh)
- gatk flagstat to check that the number of read pairs in your new file matches the old one (a bit overkill) [flagstat_metrics.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/flagstat_metrics.sh)
- get some other metrics (percent duplicates, number of reads/fraction mapped): [summary_metrics.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/summary_metrics.sh)

----
From there, the reliability of the scripts decreases. 
- unfinished script to actually call variants: [call_variants.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/call_variants.sh)
- I believe [validate_reads.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/validate_reads.sh) is a duplicate of flagstat_metrics.sh
- orphaned scripts I used simply because I was recording fixes to failed slurm jobs: [greperrors.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/greperrors.sh), [rerun_sort.sh](https://github.com/paulagardner/fox_repo/blob/data_formatting/scripts/rerun_sort.sh)
