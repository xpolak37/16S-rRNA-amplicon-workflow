process MERGING_READS {
    publishDir "${params.outdir}/bbmap", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}-mergedpairs.fastq.gz"), emit: reads
    
    script:
    """
    bbmerge.sh \\
    in1=${read1} \\
    in2=${read2} \\
    qtrim=r \\
    trimq=${params.bbmap_trimq} \
    maxlength=${params.bbmap_maxlength} \
    mininsert=${params.bbmap_mininsert} \
    threads=${task.cpus} \
    out=${sample_id}-mergedpairs.fastq.gz \
    outu=/dev/null
    """
}


process ORIENTING_READS {
    publishDir "${params.outdir}/oriented", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    path(orienting_db)
    
    output:
    tuple val(sample_id), path("${sample_id}-oriented.fastq.gz"), emit: reads
    
    script:
    """
    vsearch -orient ${reads} \\
        -db ${orienting_db} \\
        -fastqout ${sample_id}-oriented.fq \\
        -threads ${task.cpus}

    gzip -c ${sample_id}-oriented.fq > ${sample_id}-oriented.fastq.gz
    """
}