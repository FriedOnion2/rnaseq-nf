#!/usr/bin/env Rscript
# =============================================================================
#  DESeq2 差异表达分析脚本 / Differential expression analysis
#  输入：featureCounts 各样本 *counts.txt
#  输出：DEG 表格 (CSV) + 火山图 / MA 图 (PDF)
#
#  Usage: Rscript deseq2.R --outdir . --contrast treat,control
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# ---- 命令行参数（用 base R 解析，避免依赖 optparse）----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  idx <- match(paste0("--", name), args)
  if (is.na(idx)) return(default)
  args[idx + 1]
}
opts <- list(
  outdir   = get_arg("outdir", "."),
  contrast = get_arg("contrast", NULL),
  counts   = get_arg("counts", ".")
)

# ---- 读取所有样本的 raw counts ----
count_files <- list.files(opts$counts, pattern = "_counts.txt$", full.names = TRUE)
if (length(count_files) == 0) stop("未找到 *counts.txt 文件")

read_counts <- function(f) {
  sample_name <- sub("_counts.txt$", "", basename(f))
  df <- read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
                   comment.char = "#")
  # featureCounts 输出：第一列 Geneid，最后一列是该样本 count
  out <- df[, c(1, ncol(df))]
  colnames(out) <- c("Geneid", sample_name)
  out
}

counts_list <- lapply(count_files, read_counts)
cts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), counts_list)
cts[is.na(cts)] <- 0
rownames(cts) <- cts$Geneid
cts$Geneid <- NULL
cts <- round(as.matrix(cts))

# ---- 从样本名推断 condition（样本名形如 <condition>_rep<N>）----
col_data <- data.frame(
  sample    = colnames(cts),
  condition = sub("_rep\\d+$", "", colnames(cts)),
  row.names = colnames(cts),
  stringsAsFactors = FALSE
)

# 解析对比：--contrast "A,B" 表示 A(处理) vs B(对照)
# 将 condition 因子水平设为 [B, A]，使 reference = B，结果即为 A/B 的 log2FC
groups <- c("", "")
if (!is.null(opts$contrast)) {
  gr <- strsplit(opts$contrast, ",")[[1]]
  gr <- trimws(gr)
  if (length(gr) == 2) groups <- gr
}
if (nchar(groups[1]) > 0) {
  # 保证 levels 顺序为 [对照, 处理]
  lv <- unique(c(groups[2], groups[1], setdiff(unique(col_data$condition), groups)))
  col_data$condition <- factor(col_data$condition, levels = lv)
} else {
  col_data$condition <- factor(col_data$condition)
}

cat("实验设计 / design matrix:\n"); print(col_data)

# ---- 构建 DESeqDataSet ----
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design = ~ condition)

# 过滤低表达基因（至少 2 个样本 count >= 10）
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# ---- 差异表达分析 ----
dds <- DESeq(dds)

# 显式指定对比：处理组 vs 对照组
contrast_vec <- c("condition", groups[1], groups[2])
if (nchar(groups[1]) > 0) {
  cat("对比 contrast:", groups[1], "vs", groups[2], "\n")
  res <- results(dds, contrast = contrast_vec, alpha = 0.05)
} else {
  res <- results(dds, alpha = 0.05)
}

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

write.csv(res_df, file.path(opts$outdir, "deseq2_results.csv"), row.names = FALSE)

# 显著基因
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, ]
write.csv(sig, file.path(opts$outdir, "deseq2_significant.csv"), row.names = FALSE)
cat(sprintf("显著差异基因数 significant DEGs: %d\n", nrow(sig)))

# ---- 可视化 ----
# 火山图 / volcano plot
res_df$significant <- "NS"
res_df$significant[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "Up"
res_df$significant[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down"
res_df$significant <- factor(res_df$significant, levels = c("Up", "Down", "NS"))
res_df$neglog10padj <- -log10(res_df$padj)

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj, color = significant)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)",
       title = "Volcano plot") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file.path(opts$outdir, "volcano.pdf"), p_volcano, width = 7, height = 6)

# MA 图 / MA plot
p_ma <- ggplot(res_df[!is.na(res_df$padj), ], aes(x = log10(baseMean), y = log2FoldChange, color = significant)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey75")) +
  labs(x = "log10(mean normalized counts)", y = "log2 Fold Change", title = "MA plot") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file.path(opts$outdir, "ma_plot.pdf"), p_ma, width = 7, height = 6)

# 显著基因热图 / heatmap of top DEGs
# 使用 base R 的 heatmap()（无需 pheatmap 包），返回前 50 个显著 DEG 的样本×基因热图
if (nrow(sig) > 0) {
  top <- head(sig, 50)
  # 数据量小（基因数 < 1000）时 vst() 会因 nsub 默认值报错，改用 VST 直接变换
  if (nrow(dds) < 1000) {
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  } else {
    vsd <- vst(dds, blind = FALSE)
  }
  mat <- assay(vsd)[top$gene_id, , drop = FALSE]
  mat <- mat - rowMeans(mat)   # 中心化，便于看各基因相对表达
  pal <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  pdf(file.path(opts$outdir, "heatmap_top50.pdf"), width = 8, height = 10)
  heatmap(mat, col = pal, scale = "none", margins = c(10, 10),
          main = "Top 50 significant DEGs")
  dev.off()
}

cat("DESeq2 分析完成。输出文件位于:", normalizePath(opts$outdir), "\n")