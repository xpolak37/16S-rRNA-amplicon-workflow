process HOST_REMOVAL {
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}*.clean_1.fastq.gz"),
          path("${sample_id}*.clean_2.fastq.gz"), emit: reads
    
    script:
    """
    export XDG_DATA_HOME=\$PWD/.xdg
    export XDG_CACHE_HOME=\$PWD/.cache
    export NUMBA_CACHE_DIR=\${PWD}/numba_cache
    mkdir -p \${PWD}/numba_cache
    export TMPDIR=\${PWD}/tmp
    mkdir -p \${PWD}/tmp
    
    # HUMAN REMOVAL
    hostile clean \\
    --fastq1 ${read1} \\
    --fastq2 ${read2} \\
    --index human-t2t-hla-argos985-mycob140 \\
    --output . \\
    --threads ${task.cpus}


    """
}

process PHIX_REMOVAL {
    publishDir "${params.outdir}/hostile", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    path(path_bowtie_phix)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}*.clean_1.fastq.gz"),
          path("${sample_id}*.clean_2.fastq.gz"), emit: reads
    
    script:
    """
    export XDG_DATA_HOME=\$PWD/.xdg
    export XDG_CACHE_HOME=\$PWD/.cache
    export NUMBA_CACHE_DIR=\${PWD}/numba_cache
    mkdir -p \${PWD}/numba_cache
    export TMPDIR=\${PWD}/tmp
    mkdir -p \${PWD}/tmp
        
    # PHIX REMOVAL
    hostile clean \\
    --fastq1 ${read1}\\
    --fastq2 ${read2} \\
    --index ${path_bowtie_phix}/phiX174 \\
    --output . \\
    --threads ${task.cpus}
    """
}

