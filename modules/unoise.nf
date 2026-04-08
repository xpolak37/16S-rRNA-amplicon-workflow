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

    # ── 7. Map all reads back to ZOTUs to build an abundance table ──
    vsearch --usearch_global all_reads.fastq \
        --db ASV_sequences.fasta \
        --id ${params.unoise_id} \
        --fastq_qmax ${params.unoise_fastq_qmax} \
        --otutabout otu_table_raw.tsv \
        --threads ${task.cpus}

    # ── 8. Reformat: replace OTU IDs with actual sequences ──
    python3 - <<'PYEOF'
seqs = {}
with open("ASV_sequences.fasta") as f:
    header = None
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
        else:
            seqs[header] = line

with open("otu_table_raw.tsv") as f:
    lines = f.readlines()

with open("asv_table.tsv", "w") as out:
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\\t")
        if fields[0] == "#OTU ID":
            fields[0] = "SeqID"
        elif fields[0] in seqs:
            fields[0] = seqs[fields[0]]
        out.write("\\t".join(fields) + "\\n")
PYEOF
    """
}
