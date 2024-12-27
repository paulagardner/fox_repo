

params.rg_config_path = '/gpfs/data/bergstrom/foxseq2024/read-group-config.txt'
params.sequence_files_path = '/gpfs/data/bergstrom/foxseq2024'
params.index_path = '/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa'
params.output_dir = '/gpfs/data/bergstrom/paula/fox_repo/our_genomes'

//params.paired_end_fastqs = "gpfs/data/bergstrom/foxseq2024/${sample_prefix}_{1,2}.fastq.gz"


process bwa_mem_align {
    tag { "${sample_prefix}" }
    publishDir "${params.output_dir}/bwa_align", mode: 'symlink', overwrite: true
    cache 'lenient'

    input: //these are the variables that the process will actually use
        tuple val(sample_prefix), val(read_group_id), val(sample_id), val(library), val(platform), val(bamfile_basename), path(fastq1_path), path(fastq2_path)

    output:
        //path "output_dir/{bamfile_basename}.bwa.bam"  // Output to the work directory
        //path "${bamfile_basename}.bwa.bam"
	
	//set the sample ID as an output so you can carry it through the pipeline: 
        tuple val(sample_id), path("testing_${bamfile_basename}.txt"), val(bamfile_basename)
    
    script:
    """
    # Load the required modules

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

    publishDir "${params.output_dir}/sorted_files", mode: 'symlink', overwrite: true

    input:
        tuple val(sample_ID), path(bamfile), val(bamfile_basename)

    output:
        tuple val(sample_ID), path("sort_standin_${bamfile_basename}.txt"), val(bamfile_basename)

    script:
    """
    #echo "${bamfile}"
    cat "${bamfile}" > sort_standin_"${bamfile_basename}".txt
    """
}

process merge_samples {
    tag { "merge ${key}" }

    publishDir "${params.output_dir}/sorted_files", mode: 'symlink', overwrite: true

    input:
    	tuple val(key), path(list1), val(list2)

    script:
    """
    echo "Key: ${key}"
    echo "List 1: ${list1.join(', ')}"
    echo "List 2: ${list2.join(', ')}"
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
        .take(3)  // Take only the first tuple. remove these two lines, they're my equivalent of break for groovy right now.
        //.view()
        .set{first_readgroup_config_channel } // Save as a reusable channel


    //first_readgroup_config_channel.view()

    // Connect channels to processes
    bwa_output_ch = bwa_mem_align(first_readgroup_config_channel) //make bwa_output channel from bwa_mem_align
	//.view()

    sort_output_ch = sort(bwa_output_ch) //make sort_output channel
        .groupTuple()
	.view()

    sort_output_ch | merge_samples

   // merge_output = merge_samples(sort_output_ch)


   //bwa_output_ch.join(sort_output_ch)
	//.groupTuple() 
	//.view()
}


 
