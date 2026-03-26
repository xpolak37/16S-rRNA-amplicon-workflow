process HOST_REMOVAL {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)
    path(host_index_dir)

    output:
    tuple val(sample_id),
          path("*.clean_1.fastq.gz"),
          path("*.clean_2.fastq.gz"), emit: reads

    script:
    """
    export TMPDIR=\${PWD}/tmp
    mkdir -p \${PWD}/tmp
    export HOSTILE_CACHE_DIR=${host_index_dir}

    hostile clean \\
        --fastq1 ${read1} \\
        --fastq2 ${read2} \\
        --index human-t2t-hla-argos985-mycob140 \\
        --aligner bowtie2 \\
        --output . \\
        --threads ${task.cpus}
    """
}

process PHIX_REMOVAL {
    tag "$sample_id"
    publishDir "${params.outdir}/hostile", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path(path_bowtie_phix)

    output:
    tuple val(sample_id),
          path("*.clean_1.fastq.gz"),
          path("*.clean_2.fastq.gz"), emit: reads

    script:
    """
    export TMPDIR=\${PWD}/tmp
    mkdir -p \${PWD}/tmp
    export HOSTILE_CACHE_DIR=${path_bowtie_phix}

    hostile clean \\
        --fastq1 ${read1} \\
        --fastq2 ${read2} \\
        --index phiX174 \\
        --aligner bowtie2 \\
        --output . \\
        --threads ${task.cpus}
    """
}
