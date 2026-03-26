process HOST_REMOVAL {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)
    path(host_index_dir)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.host_removed.fastq.gz"),
          path("${sample_id}_R2.host_removed.fastq.gz"), emit: reads

    script:
    """
    bowtie2 \\
        -x ${host_index_dir}/human-t2t-hla-argos985-mycob140 \\
        -1 ${read1} \\
        -2 ${read2} \\
        --un-conc-gz ${sample_id}_R%.host_removed.fastq.gz \\
        --threads ${task.cpus} \\
        > /dev/null 2> bowtie2_host.log

    echo "HOST_REMOVAL summary:"
    tail -1 bowtie2_host.log
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
          path("${sample_id}_R1.decontam.fastq.gz"),
          path("${sample_id}_R2.decontam.fastq.gz"), emit: reads

    script:
    """
    bowtie2 \\
        -x ${path_bowtie_phix}/phiX174 \\
        -1 ${read1} \\
        -2 ${read2} \\
        --un-conc-gz ${sample_id}_R%.decontam.fastq.gz \\
        --threads ${task.cpus} \\
        > /dev/null 2> bowtie2_phix.log

    echo "PHIX_REMOVAL summary:"
    tail -1 bowtie2_phix.log
    """
}
