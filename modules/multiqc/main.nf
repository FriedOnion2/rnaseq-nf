// =============================================================================
//  MultiQC 模块 — 汇总所有质控报告 / Aggregate all QC reports into one HTML
// =============================================================================

process MULTIQC {
    label 'process_low'

    container 'multiqc/multiqc:v1.25.1'
    // MultiQC 镜像默认以 UID 1000 非 root 运行，无法写 root 所有的 work 目录
    containerOptions '--user 0:0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path fastqc_zips
    path trim_reports
    path star_logs
    path qualimap_reports

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/**",    emit: data, optional: true

    script:
    """
    multiqc .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}