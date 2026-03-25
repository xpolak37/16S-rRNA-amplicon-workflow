process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_${stage}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    val stage  // 'raw' or 'trimmed'
    
    output:
    path "*.html", emit: html
    path "*.zip", emit: zip
    
    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        ${params.fastqc_extra_args} \\
        ${read1} ${read2}

    """

}