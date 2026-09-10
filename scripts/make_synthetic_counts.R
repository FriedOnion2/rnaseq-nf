#!/usr/bin/env Rscript
# =============================================================================
#  生成模拟 featureCounts 输入 + 已知差异基因，用于本地验证 DESeq2 步骤。
#  Generate synthetic featureCounts-style count files with KNOWN DE genes,
#  so we can verify the DESeq2 step recovers them (ground-truth check).
#
#  Usage: Rscript scripts/make_synthetic_counts.R <outdir>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "test/counts"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
n_genes <- 1000

samples <- c("control_rep1","control_rep2","control_rep3",
             "treat_rep1","treat_rep2","treat_rep3")
group <- rep(c("control","treat"), each = 3)

# 基础表达量（对数正态分布）
base <- rlnorm(n_genes, meanlog = 6, sdlog = 0.8)  # ~ e^6 平均
base[base < 1] <- 1

# 已知差异基因：前 100 个基因在 treat 组上调 4 倍（log2FC=2），后 100 个下调
gene_id <- paste0("ENSG", sprintf("%08d", 1:n_genes))
fc <- rep(1, n_genes)
fc[1:100]   <- 4    # 上调 UP
fc[901:1000] <- 1/4 # 下调 DOWN

# 生成每样本 count 矩阵（负二项分布，DESeq2 同款假设）
# 注意：用显式向量运算而非 ifelse()，避免 ifelse 对标量条件与长向量的错误回收
for (i in seq_along(samples)) {
  disp <- 0.05                      # 较低离散度，便于干净地恢复 4x 信号
  mu <- base                        # 对照：基础表达量
  if (group[i] == "treat") {
    mu <- base * fc                 # 处理：乘上倍数变化
  }
  # 负二项抽样（size = 1/disp）
  counts <- rnbinom(n_genes, mu = mu, size = 1/disp)
  counts[counts < 0] <- 0

  df <- data.frame(
    Geneid = gene_id,
    Chr = "chr21", Start = 1, End = 100, Strand = "+",
    Length = 100,
    count = counts,
    stringsAsFactors = FALSE
  )
  # featureCounts 风格：最后一列就是该样本的 count
  colnames(df)[ncol(df)] <- samples[i]

  # 加上 featureCounts 的注释头部（真实文件里有这些）
  write.table(
    df, file.path(outdir, paste0(samples[i], "_counts.txt")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# 记录 ground truth
truth <- data.frame(gene_id = gene_id, log2FC = log2(fc), stringsAsFactors = FALSE)
write.csv(truth, file.path(outdir, "ground_truth.csv"), row.names = FALSE)

cat("生成完毕。候选差异基因真值：\n")
cat("  显著上调 (log2FC=2):  ENSG00000001 .. ENSG00000100  (100 genes)\n")
cat("  显著下调 (log2FC=-2): ENSG00000901 .. ENSG00001000 (100 genes)\n")
cat("  其余 800 个基因无差异\n")
cat("输出目录:", normalizePath(outdir), "\n")