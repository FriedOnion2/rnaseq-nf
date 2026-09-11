// =============================================================================
//  DESeq2 模块 — 差异表达分析 / Differential expression analysis (R/DESeq2)
//
//  使用公开的 Bioconductor DESeq2 镜像；本地 deseq2.R 脚本被自动整理到
//  工作目录并在容器内运行。输入为各样本 featureCounts 输出；脚本内部自动
//  合并为 raw-count 矩阵并按样本名前缀（<condition>_rep<N>）解析分组。
//
//  无需自行构建镜像。
// =============================================================================

process DESEQ2 {
    tag "deseq2"
    label 'process_low'

    container 'quay.io/biocontainers/bioconductor-deseq2:1.46.0--r44he5774e6_1'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path counts_files   // 各样本 *counts.txt（featureCounts 输出）
    val  contrast       // "A,B"（A vs B）
    path rscript        // bin/deseq2.R

    output:
    path "deseq2_results.csv",     emit: results
    path "deseq2_significant.csv", emit: significant
    path "volcano.pdf",            emit: volcano
    path "ma_plot.pdf",            emit: ma_plot
    path "heatmap_top50.pdf",      emit: heatmap
    path "*.Rout",                 emit: log, optional: true

    script:
    """
    Rscript ${rscript} \\
        --outdir . \\
        --counts . \\
        --contrast "${contrast}" 2>&1 | tee deseq2.Rout
    """

    stub:
    """
    touch deseq2_results.csv deseq2_significant.csv volcano.pdf ma_plot.pdf heatmap_top50.pdf
    """
}