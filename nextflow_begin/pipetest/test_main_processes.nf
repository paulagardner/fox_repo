

params.rg_config_path = '/gpfs/data/bergstrom/foxseq2024/read-group-config.txt'
params.sequence_files_path = '/gpfs/data/bergstrom/foxseq2024'
params.index_path = '/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa'
params.output_dir = '/gpfs/data/bergstrom/paula/fox_repo/our_genomes'

//params.paired_end_fastqs = "gpfs/data/bergstrom/foxseq2024/${sample_prefix}_{1,2}.fastq.gz"


process bwa_mem_align {
    tag { "${sample_prefix}" }
    publishDir "${params.output_dir}/bwa_align", mode: 'symlink', overwrite: true
    cache 'lenient'

    executor 'local'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue

    input: //these are the variables that the process will actually use
        tuple val(sample_prefix), val(read_group_id), val(sample_id), val(library), val(platform), val(bamfile_basename), path(fastq1_path), path(fastq2_path)

    output:
        //path "output_dir/{bamfile_basename}.bwa.bam"  // Output to the work directory
        //path "${bamfile_basename}.bwa.bam"
	
	//set the sample ID as an output so you can carry it through the pipeline: 
        tuple val(sample_id), path("testing_${bamfile_basename}.txt"), val(bamfile_basename)
    
    script:
    """

    # Ensure the output directory exists
    mkdir -p ${params.output_dir}/bwa_align 

    # Debugging output to check paths
    echo "Sample prefix: ${sample_prefix}\n \
    BAM basename: ${bamfile_basename}\n \
    Fastq 1 Path: ${fastq1_path}\n \
    Fastq 2 Path: ${fastq2_path}\n \
    Reference Index Path: ${params.index_path}\n \
    Output path: ${params.output_dir}/${bamfile_basename}.bwa.bam" \
    > 'testing_${bamfile_basename}.txt' 
    """  
}

process sort {
    tag { "Sort ${bamfile.baseName}" }
    cache 'lenient'


    executor 'local'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue

    publishDir "${params.output_dir}/sorted_files", mode: 'symlink', overwrite: true

    input:
        tuple val(sample_ID), path(bamfile), val(bamfile_basename)

    output:
        tuple val(sample_ID), path("sort_standin_${bamfile_basename}.txt"), val(bamfile_basename)

    script:
    """
    #echo "${bamfile}"
    cat "${bamfile}" > sort_standin_"${bamfile_basename}".txt
    echo "hi"
    """
}

process merge_samples {
    tag { "merge ${sample_ID}" }
    cache 'lenient'
    

    executor 'local'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue

    publishDir "${params.output_dir}/sorted_files", mode: 'symlink', overwrite: true

    input:
    	tuple val(sample_ID), path(bamfiles), val(bamfile_basenames)

    output: 
        tuple val(sample_ID), path("merged_${sample_ID}.txt"), val(bamfile_basenames)

    script:
    """
   
    #echo "Key: ${sample_ID}"
    #echo "raw list 1: ${bamfiles}"
    #echo "List 1: ${bamfiles.join(', ')}" #join is only necessary if you want to include this delimiter
    #echo "List 2: ${bamfile_basenames.join(', ')}"
    
    #cat 

    #hi
    # Concatenate all files in list1 into one file
    cat ${bamfiles.join(' ')} > merged_${sample_ID}.txt

    # Print the contents of the merged file for verification
    #echo "Merged contents of ${sample_ID}:"
    #cat merged_${sample_ID}.txt
    """
}


workflow {
    // Define the input channel
    readgroups_config_channel = Channel
        .fromPath(params.rg_config_path)
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> 
            def (sample_prefix, read_group_id, sample_id, library, platform) = row
            def bamfile_basename = sample_prefix.tokenize('/')[-1]
            def fastq1_path = "${params.sequence_files_path}/${sample_prefix}_1.fq.gz"
            def fastq2_path = "${params.sequence_files_path}/${sample_prefix}_2.fq.gz"

            return tuple(
                sample_prefix,
                read_group_id,
                sample_id,
                library,
                platform,
                bamfile_basename,
                file(fastq1_path),
                file(fastq2_path))
            }
        .set{first_readgroup_config_channel}

    // bwa_mem_align output
    bwa_output_ch = bwa_mem_align(first_readgroup_config_channel)

    sort_output_ch = sort(bwa_output_ch) //make sort_output channel
	    .groupTuple()
        .view()


    //////


    merge_output_ch = merge_samples(sort_output_ch)

}

/* workflow {
    // Define the input channel
    readgroups_config_channel = Channel
        .fromPath(params.rg_config_path)
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> 
            def (sample_prefix, read_group_id, sample_id, library, platform) = row
            def bamfile_basename = sample_prefix.tokenize('/')[-1]
            def fastq1_path = "${params.sequence_files_path}/${sample_prefix}_1.fq.gz"
            def fastq2_path = "${params.sequence_files_path}/${sample_prefix}_2.fq.gz"

            return tuple(
                sample_prefix,
                read_group_id,
                sample_id,
                library,
                platform,
                bamfile_basename,
                file(fastq1_path),
                file(fastq2_path))
            }
        .set{first_readgroup_config_channel}

    // bwa_mem_align output
    bwa_output_ch = bwa_mem_align(first_readgroup_config_channel)

    // Sort output and group by bamfile_basename
    sort_output_ch = sort(bwa_output_ch)
        .groupTuple()
        .map { tuple ->
            def (bamfile_basename, file_list, sample_list) = tuple
            def expected_file_count = sample_list.size()
            println "Checking files for ${bamfile_basename}: ${file_list.size()} files available (expected ${expected_file_count})"

            if (file_list.size() == expected_file_count) {
                // All expected files are present
                return tuple // Proceed to the merge step immediately for this group
            } else {
                // Not all files are available yet, return null to block this group
                println "Not all files available for ${bamfile_basename}. Waiting..."
                return null
            }
        }
        .filter { it != null } // Remove incomplete groups
        .view() // View for debugging

    // Create a new channel for completed groups
    completed_groups_ch = sort_output_ch

    // Merge step that processes only completed groups as they become available
    merge_output_ch = merge_samples(completed_groups_ch)
} */



/* workflow {
    // Define the input channel
    readgroups_config_channel = Channel
        .fromPath(params.rg_config_path)
        .splitCsv(skip: 1, sep: '\t')
        .map { row -> 
            def (sample_prefix, read_group_id, sample_id, library, platform) = row
            def bamfile_basename = sample_prefix.tokenize('/')[-1]
            def fastq1_path = "${params.sequence_files_path}/${sample_prefix}_1.fq.gz"
            def fastq2_path = "${params.sequence_files_path}/${sample_prefix}_2.fq.gz"

            return tuple(
                sample_prefix,
                read_group_id,
                sample_id,
                library,
                platform,
                bamfile_basename,
                file(fastq1_path),
                file(fastq2_path))
            }
        .set{first_readgroup_config_channel}

    // bwa_mem_align output
    bwa_output_ch = bwa_mem_align(first_readgroup_config_channel)

    // Sort output and group by bamfile_basename
    sort_output_ch = sort(bwa_output_ch)
        .groupTuple()
        .map { tuple ->
            def (bamfile_basename, file_list, sample_list) = tuple
            def expected_file_count = sample_list.size()
            println "Checking files for ${bamfile_basename}: ${file_list.size()} files available (expected ${expected_file_count})"

            if (file_list.size() == expected_file_count) {
                // All expected files are present
                println "Group ready for merge: ${bamfile_basename}" // Log when a group is ready
                return tuple // Proceed to the merge step immediately for this group
            } else {
                // Not all files are available yet, return null to block this group
                println "Not all files available for ${bamfile_basename}. Waiting..."
                return null
            }
        }
        .filter { it != null } // Remove incomplete groups
        //.view() // View for debugging

    // Create a new channel for completed groups
    completed_groups_ch = sort_output_ch

    // Log the channel content before passing to merge
    completed_groups_ch.view { tuple -> 
        println "Passing to merge: ${tuple[0]}" // Log which group is being passed to merge
    }

    // Merge step that processes only completed groups as they become available
    merge_output_ch = merge_samples(completed_groups_ch)
} */

//https://github.com/nextflow-io/nextflow/issues/796 //this seems to suggest that there isn't a 
//good way to \\



