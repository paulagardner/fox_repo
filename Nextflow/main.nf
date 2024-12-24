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
- The directory must contain a subdirectory named raw_sequences located two levels down.  
- A manually created readgroups file is required (automation can be explored later).
*/


/*

*/

params.rg_config_path = '/gpfs/data/bergstrom/foxseq2024/read-group-config.txt'
params.sequence_files_path = '/gpfs/data/bergstrom/foxseq2024'
params.index_path = '/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa'
params.output_dir = '/gpfs/data/bergstrom/paula/fox_repo/our_genomes'

//params.paired_end_fastqs = "gpfs/data/bergstrom/foxseq2024/${sample_prefix}_{1,2}.fastq.gz"



process sayHello {
    script:
    """
    echo 'Hello world!' 
    """
}

process test {

    publishDir "${params.output_dir}/bwa_align", mode: 'symlink'

    input:
        tuple val(sample_prefix), val(read_group_id), val(sample_id), val(library), val(platform), val(bamfile_basename), path(fastq1_path), path(fastq2_path)

    output:
        path "${bamfile_basename}.txt"

    script:
    """
    mkdir -p ${params.output_dir}/bwa_align

    # Debugging output to check paths
    echo "Sample prefix: ${sample_prefix}\n \
    BAM basename: ${bamfile_basename}\n \
    Fastq 1 Path: ${fastq1_path}\n \
    Fastq 2 Path: ${fastq2_path}\n \
    Reference Index Path: ${params.index_path}\n \
    Output path: ${params.output_dir}/${bamfile_basename}.bwa.bam" \
    > '${bamfile_basename}.txt'
    """
}


process bwa_mem_align {
    tag { "${sample_prefix}" }
    //publishDir "${params.output_dir}/bwa_align", mode: 'symlink', overwrite: true
    cache 'lenient'

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '3days'                      // Request 3 days of wall time
    memory '20 GB'                 // Request 10 GB of memory
    cpus 12                         // Request 8 CPU cores

    input: //these are the variables that the process will actually use
        tuple val(sample_prefix), val(read_group_id), val(sample_id), val(library), val(platform), val(bamfile_basename), path(fastq1_path), path(fastq2_path)

    output:
        //path "output_dir/{bamfile_basename}.bwa.bam"  // Output to the work directory
        //path "${bamfile_basename}.bwa.bam"
        tuple path("${bamfile_basename}.bwa.bam"), val(bamfile_basename)
    
    script:
    """
    # Load the required modules
    module load bwa
    module load samtools

    # Ensure the output directory exists
    mkdir -p ${params.output_dir}/bwa_align 

    # Debugging output to check paths
    echo "Sample prefix: ${sample_prefix}"
    echo "BAM basename: ${bamfile_basename}"     
    echo "Fastq 1 Path: ${fastq1_path}"
    echo "Fastq 2 Path: ${fastq2_path}"
    echo "Reference Index Path: ${params.index_path}"
    echo "Output path: ${params.output_dir}/${bamfile_basename}.bwa.bam"

    # Run bwa and samtools
    bwa mem -t ${task.cpus} -T 0 -R '@RG\\tID:${read_group_id}\\tSM:${sample_id}\\tLB:${library}\\tPL:${platform}' \
    ${params.index_path} ${fastq1_path} ${fastq2_path} | \
    samtools view -b -o "${bamfile_basename}.bwa.bam"
    """  
}

process sort {
    tag { "Sort ${bamfile.baseName}" }

    //publishDir "${params.output_dir}/sorted_files", mode: 'symlink', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '3d'                      // Request 3 days of wall time
    memory '2 GB'                 // Request 2 GB of memory
    cpus 6                         // Request 6 CPU cores


    input:
        tuple path(bamfile), val(bamfile_basename)

    output:
        tuple path("${bamfile_basename}.sort"), val(bamfile_basename)

    script:
    """
    module load samtools 
    
    samtools sort "${bamfile}" -o "${bamfile_basename}.sort"
    """
}



/*
if you need to have two input channels, this can affect resume

also, -resume/cache-ing troubleshooting: https://github.com/nextflow-io/nextflow/issues/1629
https://seqera.io/blog/demystifying-nextflow-resume/
*/

workflow {
    
    sayHello()


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
        
        //.take(2)  // Take only the first tuple. remove these two lines, they're my equivalent of break for groovy right now.
        //.set { first_readgroup_config_channel }

    // Take only the first row from the channel and pass it to bwa_mem_align process
    //readgroups_config_channel
    //    .take(1)  // Take only the first tuple
    //    .set { first_readgroup_config_channel }
    //    //.into { alignment_input_channel; test_input_channel }

    // Workflow connection
    //first_readgroup_config_channel | bwa_mem_align | sort  //COMMENT OUT TO NOT RUN THIS PROCESS DURING TESTING
    //first_readgroup_config_channel | test 

    // Connect channels to processes
    bwa_output = bwa_mem_align(readgroups_config_channel)
    sort_output = sort(bwa_output) 

}