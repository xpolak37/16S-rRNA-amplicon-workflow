process VSEARCH_UNOISE3 {
    publishDir "${params.outdir}/unoise", mode: 'copy'

    input:
    path(reads)

    output:
    path("asv_table.tsv"), emit: asv_table
    path("ASV_sequences.fasta"), emit: fasta
    path("track_control.tsv"), emit: track_control

    script:
    """
    # ── 1. Concatenate and decompress all oriented reads, labeling each read
    #       with its sample name so vsearch --otutabout produces per-sample
    #       columns.  Uses ;sample=NAME; annotation (usearch/vsearch standard)
    #       to avoid dot-splitting ambiguity.  Strip -oriented suffix.
    for f in *.fastq.gz; do
        sample=\$(basename "\$f" -oriented.fastq.gz)
        gunzip -c "\$f" | awk -v s="\$sample" '
            NR%4==1 { print "@" s ";sample=" s ";" substr(\$0,2); next }
            { print }
        '
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

    # ── 6. Assign simple ASV IDs (ASV1, ASV2, …) and strip size annotations.
    #       Simple IDs ensure the FASTA headers match the OTU table row keys,
    #       which in turn match what the taxonomy classifiers emit as SeqIDs.
    awk '/^>/{print ">ASV" ++n; next}{print}' zotus.fasta \
        | sed 's/;size=[0-9]*//' > ASV_sequences.fasta

    # ── 7. Convert labeled reads to FASTA for mapping ──
    vsearch --fastx_filter all_reads.fastq \
        --fastq_qmax ${params.unoise_fastq_qmax} \
        --fastaout all_reads.fasta

    # ── 8. Map all reads back to ZOTUs to build an abundance table ──
    vsearch --usearch_global all_reads.fasta \
        --db ASV_sequences.fasta \
        --id ${params.unoise_id} \
        --otutabout otu_table_raw.tsv \
        --threads ${task.cpus}

    # ── 9. Rename header column; ASV IDs are already the row keys ──
    awk 'BEGIN{FS=OFS="\\t"}
        NR==1 { \$1="SeqID"; print; next }
        { print }
    ' otu_table_raw.tsv > asv_table.tsv

    # ── 10. Per-step read tracking ──

    # Count input reads per sample (from FASTQ headers after step 1)
    awk 'NR%4==1 {
        s = \$0
        sub(/.*sample=/, "", s)
        sub(/;.*/, "", s)
        count[s]++
    }
    END { for (s in count) print s "\\t" count[s] }' all_reads.fastq \
        | sort > counts_input.txt

    # Count filtered reads per sample (from FASTA headers after step 2)
    awk '/^>/ {
        s = \$0
        sub(/.*sample=/, "", s)
        sub(/;.*/, "", s)
        count[s]++
    }
    END { for (s in count) print s "\\t" count[s] }' filtered.fasta \
        | sort > counts_filtered.txt

    # Sum mapped reads per sample from the OTU table (after step 8)
    awk 'BEGIN{FS=OFS="\\t"}
         NR==1 { for(i=2;i<=NF;i++) name[i]=\$i; next }
         { for(i=2;i<=NF;i++) count[i]+=\$i }
         END { for(i=2;i<=NF;i++) print name[i], count[i] }' otu_table_raw.tsv \
        | sort > counts_mapped.txt

    # Aggregate counts for pooled steps
    derep_n=\$(grep -c "^>" derep.fasta)
    zotus_raw_n=\$(grep -c "^>" zotus_raw.fasta)
    zotus_n=\$(grep -c "^>" zotus.fasta)

    # Assemble track_control.tsv
    echo -e "# Aggregate: dereplicated_uniques=\${derep_n} zotus_before_chimera=\${zotus_raw_n} zotus_after_chimera=\${zotus_n}" \
        > track_control.tsv
    echo -e "SampleID\\tinput\\tfiltered\\tmapped" >> track_control.tsv

    awk 'BEGIN{FS=OFS="\\t"}
         f==0 { inp[\$1]=\$2; samples[\$1]=1; next }
         f==1 { filt[\$1]=\$2; samples[\$1]=1; next }
         f==2 { mapped[\$1]=\$2; samples[\$1]=1; next }
         END {
             for (s in samples)
                 print s, inp[s]+0, filt[s]+0, mapped[s]+0
         }' f=0 counts_input.txt f=1 counts_filtered.txt f=2 counts_mapped.txt \
        | sort >> track_control.tsv
    """
}
