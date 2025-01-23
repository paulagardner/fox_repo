params.rg_config_path = '/gpfs/data/bergstrom/foxseq2024/read-group-config.txt'
params.bamfiles_path = '/gpfs/home/xrq24scu/fox_repo/our_genomes/results/mapping/bwamem'
params.sequence_files_path = '/gpfs/data/bergstrom/foxseq2024'
params.reference_path = '/gpfs/data/bergstrom/ref/fox/mVulVul1/bwa/mVulVul1.fa'
params.output_dir = '/gpfs/data/bergstrom/paula/fox_repo/our_genomes/results'


process calculate_coverage {
    tag { "coverage ${sample_id}" }
    cache 'lenient'
    publishDir "${params.output_dir}/bamfileinfo_test", mode: 'symlink', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '1d'                      // Request n days of wall time
    memory '16 GB'                 // Request n GB of memory
    cpus 4                         // Request n CPU cores

    input:
    	tuple val(sample_id), path(file)

    output:
        tuple val(sample_id), path("WgsMetrics.${sample_id}.txt")
        
    script:
    """
    # Ensure the output directory exists
    mkdir -p "${params.output_dir}/bamfileinfo_test"

    java -jar /gpfs/software/ada/picard/2.24.1/picard.jar CollectWgsMetrics I=$file O=WgsMetrics.${sample_id}.txt  R=${params.reference_path}

    """


}

process summarize {
    tag { "summarize ${sample_id}" }
    cache 'lenient'
    publishDir "${params.output_dir}/metrics", mode: 'symlink', overwrite: true

    // SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '1d'                      // Request n days of wall time
    memory '16 GB'                 // Request n GB of memory
    cpus 4                         // Request n CPU cores

    input:
    	tuple val(sample_id), path(file)

    output:
        tuple val(sample_id), path("WgsMetrics.${sample_id}.txt")
            
    script:
    """
    echo "##########AVERAGE COVERAGE##########" 
    for sample_id in ${sample_id}; do
        cat ${params.output_dir}/metrics/WgsMetrics.${sample_id}.txt \\
        | grep -A1 "^GENOME_TERRITORY" \\
        | cut -f 2 \\
        | sed 1d \\
        | awk -v sample="\${sample_id}" '{print sample, \$0}'
    done
    """

    // the above doesn't quite work yet, as it iterates over every sample in this file, and is not
    //responding to the sample_id being passed. However, it reads the file, so stopping here for the night. 


}

workflow {
    
    bamfile_config_channel = Channel
        .fromPath("${params.bamfiles_path}/*.bam")
        .map { file -> 
            def sample_id = file.getName().tokenize('.')[0]  // Get the part before the first `.`
            tuple(sample_id, file)  // Create a tuple containing these values s
        }
        .view()

    coverage_output_ch = calculate_coverage(bamfile_config_channel)
    summary_ch = summarize(coverage_output_ch)
}
   


