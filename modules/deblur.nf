process QIIME_IMPORT {

    input:
    path(reads)

    output:
    path("imported_data.qza")

    script:
    """

    export NUMBA_CACHE_DIR=\${PWD}/numba_cache
    mkdir -p \${PWD}/numba_cache

    # make manifest
    echo -e "sample-id\tabsolute-filepath\tdirection" > qiime_manifest.tsv
    for read in *.fastq.gz; do
        sample=\$(echo \$read | sed 's/-oriented.fastq.gz//')
        echo -e "\${sample}\t\${PWD}/\${read}\tforward" >> qiime_manifest.tsv
    done


    qiime tools import \\
        --input-path qiime_manifest.tsv \\
    	--type 'SampleData[SequencesWithQuality]' \\
    	--input-format SingleEndFastqManifestPhred33V2 \\
    	--output-path imported_data.qza
    """
} 

process QIIME_DEBLUR {
    publishDir "${params.outdir}/deblur", mode: 'copy'

    input:
    path(imported_data)

    output:
    path("asv_table.tsv"), emit: asv_table
    path("ASV_sequences.fasta"), emit: fasta
    path("stats.csv"), emit: stats

    script:
    """
export NUMBA_CACHE_DIR=\${PWD}/numba_cache
mkdir -p \${PWD}/numba_cache

export TMPDIR=\${PWD}/tmp
mkdir -p \${PWD}/tmp

# DENOISING WITH DEBLUR
qiime deblur denoise-16S \\
    --i-demultiplexed-seqs ${imported_data} \\
    --p-trim-length ${params.deblur_ptrimlength} \\
    --p-min-reads ${params.deblur_minreads} \\
    --p-min-size ${params.deblur_minsize} \\
    --p-sample-stats \\
    --o-representative-sequences ASV_seqs.qza \\
    --o-table ASV_abundance.qza \\
    --o-stats stats_deblur.qza \\
    --p-jobs-to-start ${task.cpus}

# EXPORT THE RESULTS
qiime tools export \\
    --input-path ASV_abundance.qza \\
    --output-path .

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

qiime tools export \\
    --input-path ASV_seqs.qza \\
    --output-path .

qiime tools export \\
    --input-path stats_deblur.qza \\
    --output-path .

# CONVERT THE RESULTS
sed -i '2s/^#OTU ID/OTU_ID/' feature-table.tsv

python3 - <<'EOF'
seqs = {}
with open("dna-sequences.fasta") as f:
    header = None
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
        else:
            seqs[header] = line

with open("feature-table.tsv") as f:
    lines = f.readlines()

with open("asv_table.tsv", "w") as out:
    for i, line in enumerate(lines):
        if i == 0:
            continue
        fields = line.strip().split("\\t")
        if fields[0] == "OTU_ID":
            fields[0] = "SeqID"
        elif fields[0] in seqs:
            fields[0] = seqs[fields[0]]
        out.write("\\t".join(fields) + "\\n")

with open("ASV_sequences.fasta", "w") as out_fasta:
    for seq in seqs.values():
        out_fasta.write(f">{seq}\\n{seq}\\n")
EOF
    """
}


