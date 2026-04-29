/*
========================================================================================
    MOCK COMMUNITY EVALUATION MODULE
========================================================================================
    Automatically evaluates mock community samples against reference composition.
    
    Input: MetaStandard output TSV + reference mock files
    Output: Barplot (PNG), metrics (TSV), composition table (TSV)
----------------------------------------------------------------------------------------
*/

process MOCK_EVALUATION {
    publishDir "${params.outdir}/mock_evaluation", mode: 'copy'
    
    input:
    path(metastandard_table)
    path(mock_abundance)
    path(mock_taxa)
    path(mock_synonyms)
    val(workflow_id)  // e.g., "dada2PE_NaiveBayes"

    output:
    path "*_barplot.png",      emit: barplot,     optional: true
    path "*_metrics.tsv",      emit: metrics,     optional: true
    path "*_composition.tsv",  emit: composition, optional: true

    script:
    def prefix = "${workflow_id}_mock"
    """
    python3 ${projectDir}/bin/mock_evaluation.py \\
        --metastandard_table ${metastandard_table} \\
        --mock_abundance ${mock_abundance} \\
        --mock_taxa ${mock_taxa} \\
        --synonyms ${mock_synonyms} \\
        --mock_pattern "${params.mock_pattern}" \\
        --output_prefix ${prefix} \\
        --top_n ${params.mock_top_n} \\
        --run_id ${params.run_id}
    """
}
