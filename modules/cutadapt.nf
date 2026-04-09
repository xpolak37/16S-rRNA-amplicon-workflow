// Orientation-splitting cutadapt pass used exclusively by the DADA2 paired-end branch.
// Detects reverse-oriented pairs (R_primer on R1, F_primer on R2) and emits four
// files per sample: fwd_R1/R2 (forward-oriented, primers already stripped by the
// upstream CUTADAPT call) and rev_R1/R2 (reverse-oriented, primers stripped here,
// R1/R2 swapped so DADA2 sees correct F/R polarity).
process CUTADAPT_DADA2_ORIENT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}_fwd_R1.fastq.gz"), path("${sample_id}_fwd_R2.fastq.gz"),
          path("${sample_id}_rev_R1.fastq.gz"), path("${sample_id}_rev_R2.fastq.gz")

    script:
    """
    # One cutadapt pass:
    #   -g ^R_primer / -G ^F_primer  anchored to 5' of R1 / R2 respectively
    #   -o → rev_R2, -p → rev_R1    R1/R2 are intentionally swapped here so that
    #                                 rev_R1 = physical R2 (F-primer side) and
    #                                 rev_R2 = physical R1 (R-primer side), giving
    #                                 DADA2 the correct F/R polarity for this group.
    #   --untrimmed-output/paired    forward-oriented reads pass straight through
    #                                 (primers already stripped by upstream CUTADAPT).
    cutadapt \\
        --cores ${task.cpus} \\
        -g ^${params.r_primer} -G ^${params.f_primer} \\
        -o ${sample_id}_rev_R2.fastq.gz \\
        -p ${sample_id}_rev_R1.fastq.gz \\
        --untrimmed-output ${sample_id}_fwd_R1.fastq.gz \\
        --untrimmed-paired-output ${sample_id}_fwd_R2.fastq.gz \\
        ${read1} ${read2}

    # Guard: if no reverse-oriented reads were found, cutadapt may not emit the
    # rev files at all. Create empty gz stubs so Nextflow output globs don't fail.
    [ -f ${sample_id}_rev_R1.fastq.gz ] || { printf '\\x1f\\x8b\\x08\\x00\\x00\\x00\\x00\\x00' > ${sample_id}_rev_R1.fastq.gz; }
    [ -f ${sample_id}_rev_R2.fastq.gz ] || { printf '\\x1f\\x8b\\x08\\x00\\x00\\x00\\x00\\x00' > ${sample_id}_rev_R2.fastq.gz; }
    """
}

process CUTADAPT {
    tag "$sample_id"
    publishDir "${params.outdir}/cutadapt", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: reads

    script:
    """
     # Compute reverse complements
    f_rc=\$(echo "${params.f_nextera}" | tr 'ACGTacgt' 'TGCAtgca' | rev)
    r_rc=\$(echo "${params.r_nextera}" | tr 'ACGTacgt' 'TGCAtgca' | rev)

    cutadapt \\
        --cores ${task.cpus} \\
        -g ^${params.f_primer} -G ^${params.r_primer} \\
        -a ${params.f_nextera} -A ${params.r_nextera} \\
        -A \${f_rc} -a \${r_rc} \\
        -a 'A{10}' -a 'G{10}' -g 'A{10}' -g 'G{10}' \\
        -A 'A{10}' -A 'G{10}' -G 'A{10}' -G 'G{10}' \\
        -o ${sample_id}_R1.trimmed.fastq.gz \\
        -p ${sample_id}_R2.trimmed.fastq.gz \\
        ${read1} ${read2}

    """
}