#!/usr/bin/env Rscript
# =============================================================================
#  将真实公开数据集 (airway / GEO GSE52778) 转成 featureCounts 风格输入，
#  从而用流水线里那份完全相同的 bin/deseq2.R 跑差异表达。
#
#  数据集：Himes et al. 2014 (PLoS One) 气道平滑肌 cell line RNA-seq，
#  4 个未处理 (untreated) vs 4 个地塞米松 (dexamethasone) 处理样本，
#  共 64,102 个基因。真实、可引用（GEO: GSE52778）。
#
#  输出：data/airway_counts/<sample>_counts.txt（featureCounts 风格）
#
#  Usage: Rscript scripts/make_airway_counts.R <outdir>
# =============================================================================

suppressPackageStartupMessages({ library(airway); library(SummarizedExperiment) })

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "data/airway_counts"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data(airway)
se <- airway

cnt  <- assay(se)            # 原始计数矩阵 (64,102 x 8)
meta <- as.data.frame(colData(se))

cat("样本信息 / samples:\n")
print(meta[, c("Run", "cell", "dex")])

# 样本命名：<condition>_rep<N> —— 与 bin/deseq2.R 的解析逻辑一致
# airway 中 dex 列为因子 "trt"(地塞米松处理) / "untrt"(未处理)
is_trt <- as.character(meta$dex) == "trt"
meta$sample <- NA_character_
meta$sample[is_trt]  <- paste0("treated_rep", seq_len(sum(is_trt)))
meta$sample[!is_trt] <- paste0("control_rep", seq_len(sum(!is_trt)))

gene_ids <- rownames(cnt)

for (i in seq_len(ncol(cnt))) {
  df <- data.frame(
    Geneid  = gene_ids,
    Chr     = "chrNA", Start = 1, End = 1, Strand = "+", Length = 1,
    count   = cnt[, i],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  colnames(df)[ncol(df)] <- meta$sample[i]
  write.table(df, file.path(outdir, paste0(meta$sample[i], "_counts.txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("完成。样本映射：\n")
print(meta[, c("Run", "sample")])
cat("输出目录:", normalizePath(outdir), "\n")
cat("对比参数应为 --contrast treated,control\n")