process METASTANDARD_PLOTS {
    publishDir "${params.outdir}/metastandard_plots", mode: 'copy'

    input:
    path(metastandard_tsv)

    output:
    path "*_barplot.png", emit: barplot
    path "*_heatmap.png", emit: heatmap

    script:
    def prefix = metastandard_tsv.baseName
    """
    python3 ${projectDir}/bin/plot_metastandard.py \\
        --input  ${metastandard_tsv} \\
        --top_n  ${params.metastandard_top_n} \\
        --prefix ${prefix}
    """
}
