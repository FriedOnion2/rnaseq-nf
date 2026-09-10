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
  library(pheatmap)
  library(RColorBrewer)
  library(optparse)
})

# ---- 命令行参数 ----
option_list <- list(
  make_option("--outdir",   type="character", default=".",  help="输出目录"),
  make_option("--contrast", type="character", default=NULL,  help="对比: A,B（A vs B）"),
  make_option("--counts",   type="character", default=".",  help="featureCounts 输出目录或文件列表")
)
opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailingOnly=TRUE))

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
col_data$condition <- factor(col_data$condition)

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
res <- results(dds, alpha = 0.05)

# 若指定了 contrast，据此命名输出
if (!is.null(opts$contrast)) {
  groups <- strsplit(opts$contrast, ",")[[1]]
  if (length(groups) == 2) {
    cat("对比 contrast:", groups[1], "vs", groups[2], "\n")
  }
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
if (nrow(sig) > 0) {
  top <- head(sig, 50)
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)[top$gene_id, , drop = FALSE]
  mat <- mat - rowMeans(mat)
  p_heatmap <- pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE,
                         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                         show_rownames = TRUE, main = "Top 50 significant DEGs")
  ggsave(file.path(opts$outdir, "heatmap_top50.pdf"), p_heatmap, width = 8, height = 10)
}

cat("DESeq2 分析完成。输出文件位于:", normalizePath(opts$outdir), "\n")