// =============================================================================
//  MultiQC 模块 — 汇总所有质控报告 / Aggregate all QC reports into one HTML
// =============================================================================

process MULTIQC {
    label 'process_low'

    container 'multiqc/multiqc:v1.25.1'

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
    multiqc . --outdir multiqc_data --filename multiqc_report
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html
    """
}