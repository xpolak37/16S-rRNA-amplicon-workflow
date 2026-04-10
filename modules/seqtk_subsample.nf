process SEQTK_SUBSAMPLE {
    tag "$sample_id"
    publishDir "${params.outdir}/subsampled", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.subsampled.fastq.gz"),
          path("${sample_id}_R2.subsampled.fastq.gz"),
          emit: reads

    script:
    """
    seqtk sample -s100 ${read1} ${params.quick_depth} | gzip > ${sample_id}_R1.subsampled.fastq.gz
    seqtk sample -s100 ${read2} ${params.quick_depth} | gzip > ${sample_id}_R2.subsampled.fastq.gz
    """
}
