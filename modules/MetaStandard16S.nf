process METASTANDARD {
    publishDir "${params.outdir}/metastandard", mode: 'copy'

    input:
    path(asv_table)
    path(taxa_table)
    val taxa_tool


    output:
    path "*.tsv"

    script:
    
    """
    python3 ${projectDir}/bin/metastandard16S.py \\
    --asv_table  ${asv_table} \\
    --taxa_table ${taxa_table} \\
    --taxa_tool  ${taxa_tool} \\
    --level      ${params.tax_level} \\
    --run_id     ${params.run_id}
    """
}

