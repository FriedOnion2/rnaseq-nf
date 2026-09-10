// =============================================================================
//  Qualimap 模块 — 比对后质控 / Post-alignment QC (coverage, mapping quality)
// =============================================================================

process QUALIMAP {
    tag "$sample_id"
    label 'process_medium'

    container 'quay.io/biocontainers/qualimap:2.3--hdfd78af_0'

    input:
    tuple val(sample_id), path(bam), path(gtf)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap"), emit: reports
    path "*_qualimap/**", emit: report_raw

    script:
    """
    qualimap rnaseq \\
        -bam ${bam} \\
        -gtf ${gtf} \\
        -outdir ${sample_id}_qualimap \\
        -outformat HTML \\
        --java-mem-size=4G
    """
}