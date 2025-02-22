/*
How to run this script:  
Specify paths to the directory, readgroups file, and reference genome.  
Note: Singularity may struggle with symlinks, so you might need to use full file paths.  
Alternatively, create a script to run the Nextflow command from the given directory (see README for details).

Example usage:
nextflow run main.nf --inputDir /path/to/genomeA --reference /path/to/genomeA/reference.fasta
nextflow run main.nf --reference /gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa -process.echo
nextflow run main.nf --inputDir /gpfs/data/bergstrom/paula/fox_repo/NA_genomes
*/

/*
Requirements for the sample set directory (--inputDir):  
- The directory must contain a subdirectory named raw_sequences located two levels down. ***********not sure if this is still true  
- A manually created readgroups file is required (automation can be explored later).
*/


/*

*/

params.rg_config_path = '/gpfs/data/bergstrom/foxseq2024/read-group-config.txt'
params.sequence_files_path = '/gpfs/data/bergstrom/foxseq2024'
params.reference_path = '/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa'
params.output_dir = '/gpfs/data/bergstrom/paula/fox_repo/our_genomes'

//params.paired_end_fastqs = "gpfs/data/bergstrom/foxseq2024/${sample_prefix}_{1,2}.fastq.gz"

process bwa_mem_align {
    tag { "bwa align ${bamfile_basename}" }
    publishDir "${params.output_dir}/bwa_align", mode: 'symlink', overwrite: true


    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '3days'                   // Request 3 days of wall time
    memory '20 GB'                 // Request 20 GB of memory
    cpus 12                        // Request 12 CPU cores

    input:
        tuple val(sample_prefix), val(read_group_id), val(sample_id), val(library), val(platform), val(bamfile_basename), path(fastq1_path), path(fastq2_path)

    output:
        tuple val(sample_id), path("${bamfile_basename}.bwa.bam"), val(bamfile_basename)

    script:
    """
    # Load the required modules
    module load bwa
    module load samtools
    module load seqtk

    # Ensure the output directory exists
    mkdir -p "${params.output_dir}/bwa_align"

    # Run bwa and samtools
    bwa mem -t ${task.cpus} -T 0 -R '@RG\\tID:${read_group_id}\\tSM:${sample_id}\\tLB:${library}\\tPL:${platform}' \
    ${params.reference_path} ${fastq1_path} ${fastq2_path} | \
    samtools view -b -o "${bamfile_basename}.bwa.bam"
    
    """
}



process sort {
    tag { "Sort ${bamfile.baseName}" }
    publishDir "${params.output_dir}/sorted_files", mode: 'symlink', overwrite: true


    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '3d'                      // Request n days of wall time
    memory '10 GB'                 // Request n GB of memory
    cpus 6                         // Request n CPU cores


    input:
        tuple val(sample_id), path(bamfile), val(bamfile_basename)

    output:
        tuple val(sample_id), path("${bamfile_basename}.sort"), val(bamfile_basename)

    script:
    """
    # Ensure the output directory exists
    mkdir -p "${params.output_dir}/sorted_files"

    module load samtools 
    samtools sort "${bamfile}" -o "${bamfile_basename}.sort"
    """
}


process merge_samples {
    tag { "merge ${sample_id}" }
    publishDir "${params.output_dir}/merged_files", mode: 'symlink', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '3d'                      // Request n days of wall time
    memory '4 GB'                 // Request n GB of memory
    cpus 10                         // Request n CPU cores

    input:
    	tuple val(sample_id), path(bamfiles), val(bamfiles_basenames) //these become plural because of the post-processing I 
        //do in the workflow to group the Tuple.

    output:
        tuple val(sample_id), path("${sample_id}.merged"), val(bamfiles_basenames)

    script:
    """
    # Ensure the output directory exists
    mkdir -p "${params.output_dir}/merged_files"

    module load samtools
    samtools merge ${bamfiles} -o "${sample_id}.merged"
    """
}

process mark_duplicates {
    tag { "dupmark ${sample_id}" }
    publishDir "${params.output_dir}/duplicates_marked", mode: 'copy', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '3d'                      // Request n days of wall time
    memory '16 GB'                 // Request n GB of memory
    cpus 10                         // Request n CPU cores

    input:
    	tuple val(sample_id), path(merge_file), val(bamfiles_basenames)

    output:
        tuple val(sample_id), path("${sample_id}.marked.bam"), path("metrics-MarkDuplicates.${sample_id}.txt"), val(bamfiles_basenames)
        

    script:
    """
    # Ensure the output directory exists
    mkdir -p "${params.output_dir}/duplicates_marked"

    java -jar /gpfs/software/ada/picard/2.24.1/picard.jar MarkDuplicates I=${merge_file} O=${sample_id}.marked.bam M=metrics-MarkDuplicates.${sample_id}.txt

    """
}


process validate_bamfile {
    tag { "validate bam ${sample_id}" }
    publishDir "${params.output_dir}/duplicates_marked", mode: 'copy', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '1d'                      // Request n days of wall time
    memory '16 GB'                 // Request n GB of memory
    cpus 4                         // Request n CPU cores

    input:
    	tuple val(sample_id), path(mdup_file), path(mdup_metrics), val(bamfiles_basenames)

    output:
        tuple val(sample_id), path("ValidateSamfile.${sample_id}.txt"), val(bamfiles_basenames)

    script:
    """
    # Ensure the output directory exists
    mkdir -p "{params.output_dir}/duplicates_marked"

    module load samtools
    
    java -jar /gpfs/software/ada/picard/2.24.1/picard.jar ValidateSamFile I=$mdup_file O=ValidateSamfile.${sample_id}.txt M=VERBOSE

    """
}

process index_mdup_bam {
    tag { "index bam ${sample_id}" }
    publishDir "${params.output_dir}/duplicates_marked", mode: 'copy', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '1d'                      // Request n days of wall time
    memory '16 GB'                 // Request n GB of memory
    cpus 4                         // Request n CPU cores

    input:
    	tuple val(sample_id), path(mdup_file), path(mdup_metrics), val(bamfiles_basenames)

    output:
        tuple val(sample_id), path("${mdup_file}.bai"), val(bamfiles_basenames)

    script:
    """
    # Ensure the output directory exists
    mkdir -p "{params.output_dir}/duplicates_marked"

    module load samtools
    
    samtools index -@ 4 ${mdup_file}

    """
}





/*
if you need to have two input channels, this can affect resume

also, -resume/cache-ing troubleshooting: https://github.com/nextflow-io/nextflow/issues/1629
https://seqera.io/blog/demystifying-nextflow-resume/
*/

workflow {
    
   // sayHello()


    ///////////////////////////////////////////////////////////////////////////
    // FORMATTING FOR ALIGNMENT FROM READ GROUP CONFIG FILE
    readgroups_config_channel = Channel
        .fromPath(params.rg_config_path)
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> 
            // Unpack the row into named variables
            def (sample_prefix, read_group_id, sample_id, library, platform) = row

            // Extract only the filename portion (without directory structure)
            // sample_prefix = sample_prefix.tokenize('/')[-1]  // tested with this- replace [basename_bamfile] in test and bwa mem align to replicate
            bamfile_basename = sample_prefix.tokenize('/')[-1]  // Keep only the last part-1]. hopefully, tokenize will NOT split if / does not exist, and therefore match whatever path you

            //https://www.nextflow.io/docs/latest/working-with-files.html - apparently, using ** will search through subdirectories. Could be useful...
            // Define the full paths to the FASTQ files
            //def fastq1_path = "params.sequence_files/${sample_prefix}_1.fq.gz" //may have to look into whether this will break if sequences are NOT paired-end. Consider the {1,2} format in https://github.com/nextflow-io/nextflow/discussions/2923
            def fastq1_path = "${params.sequence_files_path}/${sample_prefix}_1.fq.gz"
            def fastq2_path = "${params.sequence_files_path}/${sample_prefix}_2.fq.gz"

            // Return a tuple containing all required fields
            return tuple(
                sample_prefix,  // Sample name or prefix
                read_group_id,  // Read group ID from the config
                sample_id,  // Sample ID from the config
                library,  // Library info from the config
                platform,  // Platform info from the config
                bamfile_basename,
                file(fastq1_path),  // Path to the first FASTQ file
                file(fastq2_path),  // Path to the second FASTQ file
            )
        }
        
        //.take(5)  // Take only the first tuple. remove these two lines, they're my equivalent of break for groovy right now.
        .set { test_readgroups_config_channel }

    //test_readgroups_config_channel.take(1)

	//really should look into map() or another way to 
	// have the names of the variables assigned above, 
	// then their actual values, shown on the console or in
	// slurm out for debugging. For now, they're not labeled.

    // Connect channels to processes
    bwa_output_ch = bwa_mem_align(test_readgroups_config_channel)

    sort_output_ch = sort(bwa_output_ch)
        .groupTuple()

    /*this link seems to indicate that there is no (obvious) way around 
    the merge process waiting for ALL tasks in sort to finish. return to this later?
    */

    merge_output_ch = merge_samples(sort_output_ch)

    duplicates_output_ch = mark_duplicates(merge_output_ch)

    validation_ch  = validate_bamfile(duplicates_output_ch)

    //placeholders for processes that rely only on the duplicates marking step, and which therefore can 
    //be run as different processes with placeholders
    //coverage_output = calculate_coverage(duplicates_output_ch)

    index_ch = index_mdup_bam(duplicates_output_ch)

    //merge_output_ch.view()
    //duplicates_output_ch.view()
    //index_ch.view()
}
   