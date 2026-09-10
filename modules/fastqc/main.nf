// =============================================================================
//  FastQC 模块 — 原始测序数据质控 / Raw sequencing quality control
// =============================================================================

process FASTQC {
    tag "$sample_id"
    label 'process_low'

    container 'biocontainers/fastqc:v0.12.1_cv1'

    input:
    tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id), path("*.html"), path("*.zip"), emit: zip
    path  "*_fastqc.html"                     , emit: html
    path  "*_fastqc.zip"                      , emit: zip_only

    script:
    """
    fastqc --quiet \\
        --threads ${task.cpus} \\
        ${reads1}${reads2 ? " ${reads2}" : ''}
    """

    stub:
    """
    touch ${sample_id}_fastqc.html ${sample_id}_fastqc.zip
    """
}