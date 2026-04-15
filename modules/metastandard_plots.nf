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

process METASTANDARD_COMPARE {
    publishDir "${params.outdir}/metastandard_report", mode: 'copy'

    input:
    path(all_tsvs)

    output:
    path "compare_barplot_*.png", emit: barplots
    path "compare_braycurtis_*.png", emit: heatmaps
    path "consensus_table.tsv", emit: consensus

    script:
    """
    python3 ${projectDir}/bin/compare_metastandard.py \\
        --input  ${all_tsvs} \\
        --top_n  ${params.metastandard_top_n} \\
        --threshold 0.01 \\
        --outdir .
    """
}

process METASTANDARD_DIVERSITY {
    publishDir "${params.outdir}/metastandard_report", mode: 'copy'

    input:
    path(all_tsvs)

    output:
    path "diversity_alpha.png", emit: alpha_plot
    path "diversity_alpha.tsv", emit: alpha_table
    path "diversity_pcoa.png", emit: pcoa_plot

    script:
    """
    python3 ${projectDir}/bin/diversity_metastandard.py \\
        --input  ${all_tsvs} \\
        --outdir .
    """
}

process METASTANDARD_AGREEMENT {
    publishDir "${params.outdir}/metastandard_report", mode: 'copy'

    input:
    path(all_tsvs)

    output:
    path "agreement_jaccard_*.png", emit: jaccard_heatmaps
    path "agreement_completeness.png", emit: completeness_plot
    path "agreement_completeness.tsv", emit: completeness_table
    path "agreement_summary.tsv", emit: summary

    script:
    """
    python3 ${projectDir}/bin/agreement_metastandard.py \\
        --input  ${all_tsvs} \\
        --threshold 0.01 \\
        --outdir .
    """
}

process METASTANDARD_REPORT {
    publishDir "${params.outdir}/metastandard_report", mode: 'copy'

    input:
    path(all_files)

    output:
    path "metastandard_report.html", emit: report

    script:
    """
    python3 ${projectDir}/bin/metastandard_report.py \\
        --input  ${all_files} \\
        --outdir .
    """
}
