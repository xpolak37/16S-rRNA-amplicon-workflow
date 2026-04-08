/*
========================================================================================
    CUSTOM SUMMARY MODULE
========================================================================================
    Builds a small per-run QC report from raw FastQC outputs.

    Reports:
      - Adapter types detected across samples
      - Sequencing depth distribution on raw reads (5-number summary + mean)
      - Samples with read count below `params.min_reads_threshold`
      - Overrepresented sequences (deduped) optionally annotated with the
        top hit from a remote NCBI BLAST against `nt`. The BLAST step is
        best-effort; if it fails (no internet, NCBI rate limit, etc.) the
        report is still produced and the BLAST section is marked skipped.
----------------------------------------------------------------------------------------
*/

process CUSTOM_SUMMARY {
    publishDir "${params.outdir}/custom_summary", mode: 'copy'

    input:
    path fastqc_zips

    output:
    path "custom_summary.html",   emit: html
    path "raw_custom_summary.txt", emit: txt

    script:
    def blast_flag = params.custom_summary_blast ? "--blast" : "--no-blast"
    """
    mkdir -p fastqc_input
    for f in ${fastqc_zips}; do
        cp "\$f" fastqc_input/
    done

    python3 ${projectDir}/bin/build_custom_summary.py \\
        --fastqc-dir fastqc_input \\
        --output-html custom_summary.html \\
        --output-txt raw_custom_summary.txt \\
        --threshold ${params.min_reads_threshold} \\
        --run-id ${params.run_id} \\
        ${blast_flag}
    """
}
