params.rg_config_path = '/gpfs/data/bergstrom/foxseq2024/read-group-config.txt'
params.bamfiles_path = '/gpfs/home/xrq24scu/fox_repo/our_genomes/results/mapping/bwamem'
params.sequence_files_path = '/gpfs/data/bergstrom/foxseq2024'
params.reference_path = '/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa'
params.output_dir = '/gpfs/data/bergstrom/paula/fox_repo/our_genomes/results/variantcalling_test'



/*
#!/bin/bash

# Define input and output directories
input_dir="/gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates"
mkdir -p "$input_dir/validation_txt/slurmout"
output_dir="/gpfs/data/bergstrom/paula/fox_repo/data/marked_duplicates/validation_txt"

# Loop over each .bam file in the input directory
for bam_file in "$input_dir"/*.mdup.bam; do
	# Extract the base name, stripping out both '.merged' and '.bam' extensions if present
	base_name=$(basename "$bam_file" | sed 's/_mdup//; s/\.bam$//')

 	# Define output file paths
 	output_file="${output_dir}/ValidateSamfile.${base_name}.txt"

	  # Check if output file already exists
  	if [[ -e "$output_file" ]]; then
    		echo "Output file $output_file already exists. Skipping..."
    		continue
	fi

 	# Submit a Slurm job using --wrap
 	sbatch --job-name="bamfile_validation_${base_name}" \
        	--error "$output_dir/slurmout/${base_name}.validate.e" \
        	--output "$output_dir/slurmout/${base_name}.validate.o" \
         	--time=4-0 \
         	--mem=16G \
         	--cpus-per-task=4 \
         	--partition=compute-64-512 \
         	--wrap="java -jar /gpfs/software/ada/picard/2.24.1/picard.jar ValidateSamFile I=$bam_file O=$output_file M=VERBOSE"
	#break 
done


*/

process call_variants {
    tag { "coverage ${sample_id}" }
    cache 'lenient'
    publishDir "${params.output_dir}", mode: 'symlink', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '1d'                      // Request n days of wall time
    memory '16 GB'                 // Request n GB of memory
    cpus 4                         // Request n CPU cores

    input:
    	tuple val(sample_id), path(file)

    output:
        tuple val(sample_id), path("${sample_id}.vcf")
        
    script:
    """
    module load GATK/4.6.0.0

    gatk --java-options "-Xmx4g" HaplotypeCaller \
			-R ${params.reference_path} \
			-I ${file} \
			-O ${sample_id}.vcf \
			-ERC GVCF
    """


}

workflow {
    
    bamfile_config_channel = Channel
        .fromPath("${params.bamfiles_path}/*.bam")
        .map { file -> 
            def sample_id = file.getName().tokenize('.')[0]  // Get the part before the first `.`
            tuple(sample_id, file)  // Create a tuple containing these values s
        }
        .take(1)
        .view()
        

    variant_output_ch = call_variants(bamfile_config_channel)
}
   


   