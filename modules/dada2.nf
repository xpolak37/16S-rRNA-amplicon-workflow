process DADA2_PAIRED {
    publishDir "${params.outdir}/dada2_PE", mode: 'copy'
    
    input:
    path('*')
    
    output:
    path "asv_table.tsv", emit: asv_table
    path "track_control.tsv", emit: track_control
    path "ASV_sequences.fasta", emit: fasta
    
    script:
    """
    R1_files=\$(ls *_R1*.fastq.gz | sort | tr '\\n' ' ')
    R2_files=\$(ls *_R2*.fastq.gz | sort | tr '\\n' ' ')
    
    Rscript ${projectDir}/bin/dada2_paired.R \\
        --input_R1 \${R1_files} \\
        --input_R2 \${R2_files} \\
        --nproc ${task.cpus} \\
        --truncQ ${params.paired_truncQ} \\
        --truncLen_R1 ${params.truncLen_R1} \\
        --truncLen_R2 ${params.truncLen_R2} \\
        --maxEE_R1 ${params.maxEE_R1} \\
        --maxEE_R2 ${params.maxEE_R2} \\
        --maxMismatch ${params.maxMismatch} \\
        --minOverlap ${params.minOverlap}

    """
}

process DADA2_SINGLE {
    publishDir "${params.outdir}/dada2_SE", mode: 'copy'
    
    input:
    path(reads)
    
    output:
    path "asv_table.tsv", emit: asv_table
    path "track_control.tsv", emit: track_control
    path "ASV_sequences.fasta", emit: fasta
    
    script:
    """
    Rscript ${projectDir}/bin/dada2_single.R \\
        --input ${reads} \\
        --nproc ${task.cpus} \\
        --truncQ ${params.single_truncQ} \\
        --truncLen ${params.truncLen} \\
        --maxEE ${params.maxEE}
    """
}