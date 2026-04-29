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
    fwd_R1_files=\$(ls *_fwd_R1.fastq.gz 2>/dev/null | sort | tr '\\n' ' ')
    fwd_R2_files=\$(ls *_fwd_R2.fastq.gz 2>/dev/null | sort | tr '\\n' ' ')
    rev_R1_files=\$(ls *_rev_R1.fastq.gz 2>/dev/null | sort | tr '\\n' ' ')
    rev_R2_files=\$(ls *_rev_R2.fastq.gz 2>/dev/null | sort | tr '\\n' ' ')

    Rscript ${projectDir}/bin/dada2_paired.R \\
        --input_fwd_R1 \${fwd_R1_files} \\
        --input_fwd_R2 \${fwd_R2_files} \\
        --input_rev_R1 \${rev_R1_files} \\
        --input_rev_R2 \${rev_R2_files} \\
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