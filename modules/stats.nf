// This process computes read statistics for paired-end sequencing data using the `seqkit` tool. 
// It takes a sample ID and paths to the read files as input, and outputs a statistics file named `stats.txt`. 
// The `tag` directive is used to label the process with the sample ID for easier tracking in the workflow.
process SEQKIT_STATS {
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    path('*')
    val stage  // 'raw' or 'trimmed' ..

    output:
    path "stats_${stage}.txt", emit: stats

    script:
    """
    # Compute read statistics using seqkit
    seqkit stats --all *fastq.gz > stats_${stage}.txt
    """
}
