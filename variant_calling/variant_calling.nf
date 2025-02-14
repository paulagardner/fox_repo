params.sequence_files_path = '/gpfs/data/bergstrom/foxseq2024'
params.reference_path = '/gpfs/data/bergstrom/ref/fox/mVulVul1/mVulVul1.fa'
params.output_dir = '/gpfs/data/bergstrom/paula/fox_repo/variant_calling/outputs'
params.tsv_file = '/gpfs/data/bergstrom/paula/fox_repo/variant_calling/bamfiles.tsv'


process haplotypecaller {
    tag { "haplotypecaller ${sample_id}" }
    cache 'lenient'
    publishDir "${params.output_dir}", mode: 'symlink', overwrite: true

    //SLURM directives
    executor 'slurm'               // Use SLURM as the executor
    queue 'compute-64-512'         // Specify the SLURM partition/queue
    time '4d'                     // Request n days of wall time
    memory '32 GB'                 // Request n GB of memory
    cpus 12                   // Request n CPU cores

    input:
    	tuple val(sample_id), path(file)

    output:
        tuple val(sample_id), path("${sample_id}.g.vcf")

    script:
    """
    mkdir -p ${params.output_dir}

    module load OpenJDK/jdk-20.0.2

    java -jar /gpfs/software/ada/gatk/4.6.0.0/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar HaplotypeCaller \
		-R ${params.reference_path} \
		-I ${file} \
        -O ${sample_id}.g.vcf.gz \
		-ERC GVCF
    """
        
}

workflow {
    
    bamfile_input_channel = Channel
        .fromPath(params.tsv_file)
        .splitText()
        .map { it.trim()}
        .filter {it}
        .map { file -> 
            def sample_id = file.tokenize('/')[-2]  // Get the part before the first `.`
            tuple(sample_id, file)  // Create a tuple containing these values s
        }
        .take(1)
        .view()
    
    test_output_ch=haplotypecaller(bamfile_input_channel)
        

    ///variant_output_ch = call_variants(bamfile_config_channel)
}
   


   