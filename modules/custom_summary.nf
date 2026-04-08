/*
========================================================================================
    CUSTOM SUMMARY MODULES
========================================================================================
    Builds a small per-run QC report from raw FastQC outputs.

    Reports:
      - Adapter types detected across samples
      - Sequencing depth distribution on raw reads (5-number summary + mean)
      - Samples with read count below `params.min_reads_threshold`
      - Overrepresented sequences (deduped) optionally annotated with the
        top hit from a remote NCBI BLAST against `nt`.

    Three processes are used because no single biocontainer ships both
    Python and BLAST:
      1. CUSTOM_SUMMARY_PARSE  — Python image, parses FastQC zips
      2. CUSTOM_SUMMARY_BLAST  — BLAST image, runs remote blastn
                                  (best-effort: failure is captured into a
                                  sentinel file so the pipeline keeps going)
      3. CUSTOM_SUMMARY_RENDER — Python image, emits HTML + raw txt
----------------------------------------------------------------------------------------
*/

process CUSTOM_SUMMARY_PARSE {
    publishDir "${params.outdir}/custom_summary", mode: 'copy', pattern: 'overrepresented.fasta'

    input:
    path fastqc_zips

    output:
    path "fastqc_parsed.json",   emit: json
    path "overrepresented.fasta", emit: fasta

    script:
    """
    mkdir -p fastqc_input
    for f in ${fastqc_zips}; do
        cp "\$f" fastqc_input/
    done

    python3 ${projectDir}/bin/build_custom_summary.py parse \\
        --fastqc-dir fastqc_input \\
        --out-json fastqc_parsed.json \\
        --out-fasta overrepresented.fasta
    """
}

process CUSTOM_SUMMARY_BLAST {
    publishDir "${params.outdir}/custom_summary", mode: 'copy', pattern: 'blast_hits.tsv'

    input:
    path fasta

    output:
    path "blast_hits.tsv", emit: tsv

    script:
    """
    # Best-effort remote BLAST. If the FASTA is empty, or the network
    # call fails for any reason (no internet, NCBI rate limit, timeout,
    # …), we write a sentinel file instead of failing the pipeline.
    if [ ! -s ${fasta} ]; then
        echo "BLAST_SKIPPED: no overrepresented sequences" > blast_hits.tsv
        exit 0
    fi

    set +e
    timeout ${params.custom_summary_blast_timeout} blastn \\
        -query ${fasta} \\
        -db nt \\
        -remote \\
        -outfmt '6 qseqid ssaccver stitle bitscore' \\
        -out blast_raw.tsv 2> blast.err
    rc=\$?
    set -e

    if [ \$rc -ne 0 ]; then
        echo "BLAST_SKIPPED: blastn exited with status \$rc" > blast_hits.tsv
        if [ -s blast.err ]; then
            echo "---- blastn stderr ----" >> blast_hits.tsv
            cat blast.err >> blast_hits.tsv
        fi
        exit 0
    fi

    # Keep only the best hit per query (highest bitscore), matching the
    # original script logic. Output: qseqid <tab> accession description
    awk '
    {
        qid=\$1; acc=\$2; bs=\$NF
        # stitle is everything between field 2 and the last field
        stitle=""
        for(i=3;i<NF;i++) stitle = stitle (i==3?"":OFS) \$i
        if (!(qid in best) || bs > best[qid]) {
            best[qid] = bs
            line[qid] = qid "\\t" acc " " stitle
        }
    }
    END { for (q in line) print line[q] }
    ' blast_raw.tsv > blast_hits.tsv
    """
}

process CUSTOM_SUMMARY_RENDER {
    publishDir "${params.outdir}/custom_summary", mode: 'copy'

    input:
    path parsed_json
    path blast_tsv

    output:
    path "custom_summary.html",   emit: html
    path "raw_custom_summary.txt", emit: txt

    script:
    def attempted_flag = params.custom_summary_blast ? "--blast-attempted" : "--no-blast-attempted"
    """
    python3 ${projectDir}/bin/build_custom_summary.py render \\
        --in-json ${parsed_json} \\
        --blast-tsv ${blast_tsv} \\
        --output-html custom_summary.html \\
        --output-txt raw_custom_summary.txt \\
        --threshold ${params.min_reads_threshold} \\
        --run-id ${params.run_id} \\
        ${attempted_flag}
    """
}
