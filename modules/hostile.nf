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

process FASTQ_SYNC {
    tag "$sample_id"
    publishDir "${params.outdir}/hostile", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.decontam_synced.fastq.gz"),
          path("${sample_id}_R2.decontam_synced.fastq.gz"), emit: reads
    path("${sample_id}_fastp.json"), emit: json
    path("${sample_id}_fastp.html"), emit: html

    script:
    """
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${sample_id}_R1.decontam_synced.fastq.gz \\
        --out2 ${sample_id}_R2.decontam_synced.fastq.gz \\
        --length_required 2 \\
        --disable_adapter_trimming \\
        --disable_quality_filtering \\
        --thread ${task.cpus} \\
        --json ${sample_id}_fastp.json \\
        --html ${sample_id}_fastp.html
    """
}
