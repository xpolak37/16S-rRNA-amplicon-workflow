process RENAME_SAMPLES {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(new_id), path("${r1name}"), path("${r2name}"), emit: reads

    script:
    new_id = sample_id.replace('_', '-')
    // keep the original suffix
    // swap only the sample-ID prefix's underscores for hyphens
    r1name = new_id + r1.name.substring(sample_id.length())
    r2name = new_id + r2.name.substring(sample_id.length())
    """
    cp -L ${r1} ${r1name}.tmp && mv ${r1name}.tmp ${r1name}
    cp -L ${r2} ${r2name}.tmp && mv ${r2name}.tmp ${r2name}
    """
}