// pull read group information from the fasta header as 
// well as the eager-style input file and add it to the
// read group information

// @RG\tID:foo\tLB:lib1\tPL:bar\tSM:sample\tPU:unit1 
// ID: doesn't really seem like there is an equivalent 
// in the eager output, ask anders'advice 
// LB: Library, from input tsv
// PL: Platform, from input tsv
// SM: Sample, from input tsv
// PU: Platform Unit, from fasta--- lane, in our case, I think

// to run the script, module load nextflow/25.04.6 on HALI as of 11/02/2026
// then, run <nextflow run <scriptname>.nf --inputtsv eager_samples.tsv>
// example run with one genome: nextflow run parabricks.nf --inputtsv Lyu_genomes/inputFile.tsv 


/// locally:  nextflow run anders_rg_parabricks.nf --inputtsv multiinputtest.tsv 

//// or via: 
// mkdir -p nextflow_logs &&
// sbatch \
//   --job-name=nextflow_test \
//   --output=nextflow_logs/%x.%j.out \
//   --error=nextflow_logs/%x.%j.err \
//   --time=1:00:00 \
//   --wrap "module load nextflow/25.04.6; \
//           cd /gpfs/home/xrq24scu/fox_repo/read_mapping; \
//           nextflow run anders_rg_parabricks.nf --inputtsv multiinputtest.tsv"


// make sure to transfer to an actual container you start to build  


nextflow.enable.dsl=2


// TEMPORARY TESTING SETTINGS — REMOVE FOR PRODUCTION
params.test_mode      = true
params.test_n_reads   = 10

Channel
    .fromPath(params.inputtsv)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        tuple(
            row.Sample_Name,
            row.Library_ID,
            row.Lane,
            file(row.R1),
            file(row.R2)
        )
    }
    .set { input_ch }




// process make_readgroup {

//     tag "$Sample_Name"

//     container "/gpfs/home/xrq24scu/fox_repo/containers/R_base.sif"

//     input:
//     tuple val(Sample_Name), val(Library_ID), val(Lane), path(R1), path(R2)

//     output:
//     tuple val(Sample_Name), val(Library_ID), val(Lane), path(R1), path(R2), val(readgroup)

//     script:
    
//     """
//     readgroup=\$(python3 /gpfs/data/bergstrom/sw/bin/fixrg.py \
//         --fastq $R1 \
//         --format-bwa \
//         --tag SM:${Sample_Name} \
//         --tag LB:${Library_ID} \
//         --tag PL:ILLUMINA \
//         --tag PU:${Lane} \
//         --subsample ${params.test_n_reads})

//     echo "\$readgroup"
//     """
// }

process make_readgroup {

    tag "$Sample_Name"


    //should really reacreate this container to be a .sif and not just unspecified'
    //also, this isn't actually what runs the container logic, specifying it below works
    //container "/gpfs/data/bergstrom/paula/fox_repo/containers/R_base"
    
    // executor 'slurm'         // <<< FORCE Slurm execution
    // cpus 1
    // memory '4 GB'
    // time '1h'

    input:
    tuple val(Sample_Name), val(Library_ID), val(Lane), path(R1), path(R2)

    output:
    tuple val(Sample_Name), val(Library_ID), val(Lane), path(R1), path(R2), stdout

    script:
    """
    #trying to enable 
    module load apptainer

    echo "[CONTAINER CHECK] hostname: \$(hostname)" >&2
    echo "[CONTAINER CHECK] APPTAINER_NAME: \${APPTAINER_NAME:-not_in_container}" >&2
    echo "[CONTAINER CHECK] APPTAINER_CONTAINER: \${APPTAINER_CONTAINER:-not_in_container}" >&2

    #add syntax to detect if it's an SRR file and assign readgroups another way  

    # Run fixrg.py on a subsample to get read groups
    apptainer exec /gpfs/home/xrq24scu/fox_repo/containers/R_base \
        python3 /gpfs/data/bergstrom/sw/bin/fixrg.py \
            --fastq $R1 \
            --format-bwa \
            --tag SM:${Sample_Name} \
            --tag LB:${Library_ID} \
            --tag PL:ILLUMINA \
            --tag PU:${Lane} \
            --format-parabricks-pe-list \
            --subsample ${params.test_n_reads}  
    """
}




process parabricks_fq2bam {

    tag "$Sample_Name"

    executor 'slurm'
    queue 'gpu'

    // cpus 6
    // memory '480 GB'
    // time '12h'

    cpus 1
    memory '4 GB'
    time '15m'


    clusterOptions '--qos=gpu --gpus=1'

    publishDir "results", mode: 'symlink'

    input:
    tuple val(Sample_Name), val(Library_ID), val(Lane), path(R1), path(R2), val(readgroup)

    output:
    path "${Sample_Name}.bam"

    script:
    """
    set -euo pipefail

    module load apptainer
    nvidia-smi

    mkdir -p /gpfs/scratch/xrq24scu/${Sample_Name}

    #############################################
    # Always subsample a small number of reads for a test run
    echo "Subsampling FASTQs for minimal test run"

    zcat $R1 | head -n 400 | gzip > test_R1.fastq.gz
    zcat $R2 | head -n 400 | gzip > test_R2.fastq.gz

    # Overwrite R1/R2 for the command
    R1=test_R1.fastq.gz
    R2=test_R2.fastq.gz

    # Escape for Groovy string interpolation
    RG_STRING=\$(echo "${readgroup}" | tr '\\n' ' ')


    echo "Using readgroups:"
    echo "\$RG_STRING"
    ##############################################


    apptainer exec --nv /gpfs/data/bergstrom/paula/fox_repo/read_mapping/parabricks-4.2.0-1.sif \\
      pbrun fq2bam \\
      --num-gpus 1 \\
      --num-cpu-threads 6 \\
      --tmp-dir /gpfs/scratch/xrq24scu/${Sample_Name} \\
      --logfile ${Sample_Name}.parabricks.log \\
      --ref /gpfs/home/xrq24scu/fox_repo/read_mapping/mVulVul1/mVulVul1.fa \\
      --in-fq \$R1 \$R2 "\$RG_STRING" \\
      --out-bam ${Sample_Name}.bam
    """
}








// process parabricks_fq2bam {

//     tag "$Sample_Name"

//     executor 'slurm'
//     queue 'gpu'

//     cpus 6
//     memory '480 GB'
//     time '12h'

//     clusterOptions '--qos=gpu --gpus=1'

//     publishDir "results", mode: 'symlink'

//     input:
//     tuple val(Sample_Name), val(Library_ID), val(Lane), path(R1), path(R2), val(readgroup)

    
//     output:
//     path "${Sample_Name}.bam"

//     script:
//     """
//     set -euo pipefail

//     module load apptainer
//     nvidia-smi

//     mkdir -p /gpfs/scratch/xrq24scu/${Sample_Name}
    
//     echo "Using readgroup: ${readgroup}"

//     apptainer run --nv /gpfs/data/bergstrom/paula/fox_repo/read_mapping/parabricks-4.2.0-1.sif \\
//       pbrun fq2bam \\
//       --num-gpus 1 \\
//       --num-cpu-threads 6 \\
//       --tmp-dir /gpfs/scratch/xrq24scu/${Sample_Name} \\
//       --logfile ${Sample_Name}.parabricks.log \\
//       --ref /gpfs/home/xrq24scu/fox_repo/read_mapping/mVulVul1/mVulVul1.fa \\
//       --in-fq $R1 $R2 "${readgroup.replaceAll('\n','')}" \\
//       --out-bam ${Sample_Name}.bam
//     """
// }


workflow {

    readgroup_ch = make_readgroup(input_ch)

    readgroup_ch
        .view { Sample_Name, Library_ID, Lane, R1, R2, readgroup ->
            log.info "[RG] ${Sample_Name}: ${readgroup}"
        }

    parabricks_fq2bam(readgroup_ch)
}