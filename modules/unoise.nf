process VSEARCH_UNOISE3 {
    publishDir "${params.outdir}/unoise", mode: 'copy'

    input:
    path(reads)

    output:
    path("asv_table.tsv"), emit: asv_table
    path("ASV_sequences.fasta"), emit: fasta
    path("unoise_track_control.tsv"), emit: track_control

    script:
    """
    # ── 1. Concatenate and decompress all oriented reads, labeling each read
    #       with its sample name so vsearch --otutabout produces per-sample
    #       columns (sample name = filename minus .fastq.gz suffix).
    for f in *.fastq.gz; do
        sample=\$(basename "\$f" .fastq.gz)
        gunzip -c "\$f" | awk -v s="\$sample" '
            NR%4==1 { print "@" s "." substr(\$0,2); next }
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

    # ── 10. Build per-sample read tracking table ──
    python3 - <<'PYEOF'
import os, gzip
from collections import Counter

# Count reads per sample from input fastq.gz files
input_counts = Counter()
for f in sorted(os.listdir(".")):
    if f.endswith(".fastq.gz") and f != "all_reads.fastq":
        sample = f.replace(".fastq.gz", "")
        n = 0
        with gzip.open(f, "rt") as fh:
            for line in fh:
                if line.startswith("@"):
                    n += 1
        input_counts[sample] = n

# Count reads per sample from filtered.fasta (labels: sample.original_id)
filtered_counts = Counter()
with open("filtered.fasta") as fh:
    for line in fh:
        if line.startswith(">"):
            sample = line[1:].split(".")[0]
            filtered_counts[sample] += 1

# Final read counts from ASV table (column sums)
final_counts = Counter()
with open("asv_table.tsv") as fh:
    header = fh.readline().strip().split("\\t")
    samples = header[1:]
    for line in fh:
        fields = line.strip().split("\\t")
        for s, v in zip(samples, fields[1:]):
            final_counts[s] += int(float(v))

# Count ZOTUs before and after chimera removal
zotus_raw = sum(1 for l in open("zotus_raw.fasta") if l.startswith(">"))
zotus_clean = sum(1 for l in open("zotus.fasta") if l.startswith(">"))

all_samples = sorted(input_counts.keys())
with open("unoise_track_control.tsv", "w") as out:
    out.write("SampleID\\tinput\\tfiltered\\tmapped\\tzotus_raw\\tzotus_nonchim\\n")
    for s in all_samples:
        out.write(f"{s}\\t{input_counts[s]}\\t{filtered_counts[s]}\\t{final_counts.get(s,0)}\\t{zotus_raw}\\t{zotus_clean}\\n")
PYEOF
    """
}
