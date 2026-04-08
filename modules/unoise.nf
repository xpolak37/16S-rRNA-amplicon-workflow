process VSEARCH_UNOISE3 {
    publishDir "${params.outdir}/unoise", mode: 'copy'

    input:
    path(reads)

    output:
    path("asv_table.tsv"), emit: asv_table
    path("ASV_sequences.fasta"), emit: fasta

    script:
    """
    # ── 1. Concatenate and decompress all oriented reads ──
    for f in *.fastq.gz; do
        gunzip -c "\$f"
    done > all_reads.fastq

    # ── 2. Quality filtering ──
    vsearch --fastq_filter all_reads.fastq \
        --fastq_maxee ${params.unoise_maxee} \
        --fastq_minlen ${params.unoise_minlen} \
        --fastq_qmax ${params.unoise_fastq_qmax} \
        --fastaout filtered.fasta \
        --threads ${task.cpus}

    # ── 3. Dereplicate ──
    vsearch --derep_fulllength filtered.fasta \
        --output derep.fasta \
        --sizeout \
        --minuniquesize ${params.unoise_minsize} \
        --threads ${task.cpus}

    # ── 4. Denoise with UNOISE3 ──
    vsearch --cluster_unoise derep.fasta \
        --unoise_alpha ${params.unoise_alpha} \
        --minsize ${params.unoise_minsize} \
        --centroids zotus_raw.fasta \
        --threads ${task.cpus}

    # ── 5. De novo chimera removal (uchime3) ──
    vsearch --uchime3_denovo zotus_raw.fasta \
        --nonchimeras zotus.fasta

    # ── 6. Strip size annotations from ZOTU headers ──
    sed 's/;size=[0-9]*//' zotus.fasta > ASV_sequences.fasta

    # ── 7. Convert all reads to FASTA (usearch_global does not accept --fastq_qmax) ──
    vsearch --fastx_filter all_reads.fastq \
        --fastq_qmax ${params.unoise_fastq_qmax} \
        --fastaout all_reads.fasta

    # ── 8. Map all reads back to ZOTUs to build an abundance table ──
    vsearch --usearch_global all_reads.fasta \
        --db ASV_sequences.fasta \
        --id ${params.unoise_id} \
        --otutabout otu_table_raw.tsv \
        --threads ${task.cpus}

    # ── 9. Reformat: replace ZOTU IDs with actual sequences ──
    awk 'BEGIN{FS=OFS="\\t"}
        # Pass 1: build header→sequence map from FASTA
        NR==FNR {
            if (substr(\$0,1,1)==">") { h=substr(\$0,2); next }
            seqs[h]=seqs[h] \$0
            next
        }
        # Pass 2: rewrite the OTU table
        {
            if (\$1=="#OTU ID") { \$1="SeqID" }
            else if (\$1 in seqs) { \$1=seqs[\$1] }
            print
        }
    ' ASV_sequences.fasta otu_table_raw.tsv > asv_table.tsv
    """
}
