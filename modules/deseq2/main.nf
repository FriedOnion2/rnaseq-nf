// =============================================================================
//  DESeq2 模块 — 差异表达分析 / Differential expression analysis (R/DESeq2)
//
//  使用本地构建的镜像（modules/deseq2/Dockerfile），内含 deseq2.R 脚本。
//  输入为各样本 featureCounts 输出；脚本内部自动合并为 raw-count 矩阵并按
//  样本名前缀（<condition>_rep<N>）解析分组，再执行 DESeq2 与可视化。
//
//  构建镜像: docker build -f modules/deseq2/Dockerfile -t rnaseq-nf/deseq2:1.0 .
// =============================================================================

process DESEQ2 {
    label 'process_low'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    container 'rnaseq-nf/deseq2:1.0'

    input:
    path counts_files   // 各样本 *counts.txt（featureCounts 输出，含文件头）
    val  contrast       // "A,B"（A vs B）

    output:
    path "deseq2_results.csv",     emit: results
    path "deseq2_significant.csv", emit: significant
    path "volcano.pdf",            emit: volcano
    path "ma_plot.pdf",            emit: ma_plot
    path "heatmap_top50.pdf",      emit: heatmap
    path "*.Rout",                 emit: log, optional: true

    script:
    """
    Rscript /opt/deseq2.R \\
        --outdir . \\
        --counts . \\
        --contrast ${contrast} 2>&1 | tee deseq2.Rout
    """
}