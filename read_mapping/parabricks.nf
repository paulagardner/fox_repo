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


//// or via: 
// sbatch --job-name=nextflow_test \
//        --output=nextflow_logs/%x.%j.out \
//        --error=nextflow_logs/%x.%j.err \
//        --cpus-per-task=4 \
//        --mem=16G \
//        --time=2:00:00 \
//        --wrap "module load nextflow/25.04.6; nextflow run parabricks.nf --inputtsv Lyu_genomes/inputFile.tsv -resume"


// make sure to transfer to an actual container you start to build  


// nextflow.enable.dsl=2

// Channel
//     .fromPath(params.inputtsv)
//     .splitCsv(header: true, sep: '\t')
//     .map { row ->
//         tuple(
//             row.Sample_Name,
//             row.Library_ID,
//             row.Lane,
//             row.SeqType,
//             row.R1,
//             row.R2
//         )
//     }
//     .set { input_ch }


// process make_readgroup {

//     tag "$sample"

//     input:
//     tuple val(sample), val(library), val(lane), val(seqtype), path(r1), path(r2)

//     output:
//     tuple val(sample), val(rg)

//     script:
//     rg = "@RG\tID:${sample}\tSM:${sample}\tLB:${library}\tPL:${seqtype}\tPU:${lane}"
//     """
//     echo "$rg"
//     """
// }


// workflow {

//     make_readgroup(input_ch)
//         .view { sample, rg ->
//             log.info "[RG] ${sample}: ${rg}"
//         }

//     // parabricks_fq2bam(readgroup_ch)
// }

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



nextflow.enable.dsl=2


// TEMPORARY TESTING SETTINGS — REMOVE FOR PRODUCTION
params.test_mode      = true
params.test_n_reads   = 10000



Channel
    .fromPath(params.inputtsv)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        tuple(
            row.Sample_Name,
            row.Library_ID,
            row.Lane,
            row.SeqType,
            row.R1,
            row.R2
        )
    }
    .set { input_ch }



process make_readgroup {

    tag "$Sample_Name"

    input:
    tuple val(Sample_Name), val(Library_ID), val(Lane), val(SeqType), path(R1), path(R2)

    output:
    tuple val(Sample_Name), val(readgroup), path(R1), path(R2)

    script:
    readgroup = "@RG\tID:${Sample_Name}\tSM:${Sample_Name}\tLB:${Library_ID}\tPL:NA\tPU:${Lane}"
    """
    printf '%s' "$readgroup"
    """
}

// process subset_fastqs {

//     tag "$Sample_Name"

//     // run on Slurm, small resources since testing
//     // executor 'slurm'
//     // cpus 1
//     // memory '4 GB'
//     // time '10m'

//     input:
//     tuple val(Sample_Name), val(readgroup), path(R1), path(R2)

//     output:
//     tuple val(Sample_Name), val(readgroup),
//           path("${Sample_Name}_R1.sub.fq.gz"),
//           path("${Sample_Name}_R2.sub.fq.gz")

//     script:
//     """

//     zcat SRR15858295_1.fastq.gz | head -n 40000 | gzip > ZH_genome_ngs_R1.sub.fq.gz
//     zcat SRR15858295_2.fastq.gz | head -n 40000 | gzip > ZH_genome_ngs_R2.sub.fq.gz

//     """
// }

process parabricks_fq2bam {

    tag "$Sample_Name"

    executor 'slurm'
    queue 'gpu'
    cpus 32
    memory '128 GB'
    time '12h'

    clusterOptions '--partition=gpu --qos=gpu --gpus=1'

    input:
    tuple val(Sample_Name), val(readgroup), path(R1), path(R2)

    script:
    """
    set -euo pipefail

    echo "SLURM_JOB_ID=\$SLURM_JOB_ID"
    echo "SLURM_GPUS=\$SLURM_GPUS"
    echo "CUDA_VISIBLE_DEVICES=\$CUDA_VISIBLE_DEVICES"

    module load apptainer
    nvidia-smi

    apptainer run --nv /gpfs/data/bergstrom/paula/fox_repo/read_mapping/parabricks-4.2.0-1.sif \\
        pbrun fq2bam \\
        --num-gpus 1 \\
        --ref /gpfs/home/xrq24scu/fox_repo/read_mapping/mVulVul1/mVulVul1.fa \\
        --in-fq $R1 $R2 "$readgroup" \\
        --out-bam ${Sample_Name}.bam
    """
}


workflow {

    // readgroup_ch = make_readgroup(input_ch)
    //     .view { Sample_Name, readgroup, R1, R2 ->
    //         log.info "[RG] ${Sample_Name}: ${readgroup}"
    //     }

    // parabricks_fq2bam(readgroup_ch)


    readgroup_ch = make_readgroup(input_ch)

    // test_ch = params.test_mode \
    //     ? subset_fastqs(readgroup_ch) \
    //     : readgroup_ch
    // parabricks_fq2bam(test_ch)

    parabricks_fq2bam(readgroup_ch)
} 