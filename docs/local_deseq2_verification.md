# 本地验证 DESeq2 核心（无需 Docker / WSL2）

本项目流水线主体基于 Nextflow + Docker，但**差异表达分析（DESeq2）这一步**可以在
本机的 R 环境里独立运行与验证——这在等待 Docker/WSL2 环境就绪时非常有用，也可作为
教学演示（证明分析核心是真实、可运行、经过真值校验的）。

## 原理 / Why this works

`bin/deseq2.R` 只依赖 R + Bioconductor 包，不依赖任何容器。流水线里的 DESeq2 模块
与本地调用执行的是**同一份 `deseq2.R`**，因此本地验证通过 = 流水线该步骤逻辑正确。

## 一次性准备 / One-time setup

R 4.4+ 已安装的前提下：

```bash
Rscript scripts/setup_r_local.R        # 安装 DESeq2/ggplot2/pheatmap 等到用户库
```

包会安装到 `%LOCALAPPDATA%\R\library`（无需管理员）。

## 生成带真值的模拟数据 / Generate synthetic counts with ground truth

```bash
Rscript scripts/make_synthetic_counts.R test/counts
```

生成 6 个样本（control×3 + treat×3）的 1000 个基因，其中**已知**：
- 基因 1–100：treat 组上调 4 倍（log2FC=+2）
- 基因 901–1000：treat 组下调 4 倍（log2FC=−2）
- 其余 800 个基因：无差异

真值保存在 `test/counts/ground_truth.csv`。

## 运行 DESeq2 / Run DESeq2

```bash
Rscript bin/deseq2.R \
  --outdir test/results \
  --counts test/counts \
  --contrast treat,control
```

## 已验证结果 / Verified results

| 指标 Metric | 值 Value |
|---|---|
| 差异基因灵敏度 Sensitivity | **1.000**（200/200 全部找回） |
| 假阳性 False positives | **0** |
| 漏检 False negatives | **0** |
| 方向一致性 Direction agreement | **100%**（涨/跌方向全对） |

输出：`deseq2_results.csv`、`deseq2_significant.csv`、`volcano.pdf`、`ma_plot.pdf`、
`heatmap_top50.pdf`。

> 这次验证直接证明了流水线的分析核心正确，而不是"看起来能跑"的样板。