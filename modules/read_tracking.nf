process READ_TRACKING {
    publishDir "${params.outdir}/read_tracking", mode: 'copy'

    input:
    path(all_logs)

    output:
    path "read_tracking.tsv", emit: table
    path "read_tracking.png", emit: plot

    script:
    """
    python3 ${projectDir}/bin/read_tracking.py \\
        --input  ${all_logs} \\
        --outdir .
    """
}
