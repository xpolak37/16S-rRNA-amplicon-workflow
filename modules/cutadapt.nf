process CUTADAPT {
    tag "$sample_id"
    publishDir "${params.outdir}/cutadapt", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: reads

    script:
    """
     # Compute reverse complements
    f_rc=\$(echo "${params.f_nextera}" | tr 'ACGTacgt' 'TGCAtgca' | rev)
    r_rc=\$(echo "${params.r_nextera}" | tr 'ACGTacgt' 'TGCAtgca' | rev)

    cutadapt \\
        --cores ${task.cpus} \\
        -g ^${params.f_primer} -G ^${params.r_primer} \\
        -a ${params.f_nextera} -A ${params.r_nextera} \\
        -A \${f_rc} -a \${r_rc} \\
        -a 'A{10}' -a 'G{10}' -g 'A{10}' -g 'G{10}' \\
        -A 'A{10}' -A 'G{10}' -G 'A{10}' -G 'G{10}' \\
        -o ${sample_id}_R1.trimmed.fastq.gz \\
        -p ${sample_id}_R2.trimmed.fastq.gz \\
        ${read1} ${read2}

    """
}