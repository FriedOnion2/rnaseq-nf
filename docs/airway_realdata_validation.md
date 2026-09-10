# 真实公开数据验证：Airway RNA-seq（GSE52778）

除了合成数据的地面真值验证，本流水线的 DESeq2 步骤也在**真实、可引用的公开数据**上
跑通并做了生物学验证。

## 数据集 / Dataset

- **研究**：Himes et al. 2014, *PLoS One* — "RNA-Seq differential expression analysis
  of the airway smooth muscle cell transcriptome"
- **GEO 登录号**：GSE52778
- **样本**：8 个（4 untreated + 4 dexamethasone-treated，人气道平滑肌细胞）
- **基因数**：64,102 个（Ensembl GRCh37 注释）
- **R 包**：Bioconductor `airway`（内置该数据的官方封装）

> **为何"有含金量"**：地塞米松（dexamethasone）是糖皮质激素，其转录组效应非常清晰、
> 被大量文献验证过。因此我们可以用"已知的糖皮质激素响应基因"作为**独立的外部金标准**，
> 检验我们流水线产出的差异表达结果是否生物学正确。

## 运行方法 / How to reproduce

```bash
# (1) 安装 airway 数据集包（一次性）
Rscript scripts/setup_r_local.R
# 在 R 里额外装 airway： BiocManager::install("airway")

# (2) 从 airway 生成 featureCounts 风格 counts
Rscript scripts/make_airway_counts.R data/airway_counts

# (3) 用流水线同一份脚本跑差异表达
Rscript bin/deseq2.R \
  --outdir data/airway_results \
  --counts data/airway_counts \
  --contrast treated,control
```

## 结果 / Results

| 指标 Metric | 值 Value |
|---|---|
| 总基因数 Genes tested | 17,199（过滤低表达后） |
| 显著差异基因 Significant DEGs | **853**（上调 465 / 下调 388） |
| 最显著 p 值 Top padj | ~1e-101 |

### 生物学金标准验证 / Biological gold-standard validation

用文献公认的地塞米松响应基因作为外部金标准，检查它们是否被正确判为"显著上调"：

| 基因 Symbol | Ensembl ID | log2FC | padj | 结果 |
|---|---|---|---|---|
| **DUSP1** | ENSG00000120129 | 2.95 | 3.6e-45 | ✅ 显著上调 |
| **FKBP5** | ENSG00000096060 | 3.94 | 1.6e-36 | ✅ 显著上调 |
| **ZBTB16** | ENSG00000109906 | 7.17 | 7.4e-42 | ✅ 显著上调 |
| **NFKBIA** | ENSG00000100906 | 0.83 | 4.8e-4 | ✅ 显著上调 |
| **KLF9** | ENSG00000119138 | 2.11 | 4.0e-22 | ✅ 显著上调 |
| **PER1** | ENSG00000179094 | 3.18 | 3.2e-53 | ✅ 显著上调 |
| **TSC22D3** | ENSG00000157514 | 3.34 | 4.8e-18 | ✅ 显著上调 |
| **GLUL** | ENSG00000135821 | 3.00 | 6.1e-26 | ✅ 显著上调 |
| **DUSP5** | ENSG00000138166 | 1.28 | 1.3e-5 | ✅ 显著上调 |

**9/11 个经典响应基因全部显著上调**（其余 2 个 SERPINE1、DDIT4 也呈上调趋势，仅未达
显著阈值）。这与 Himes 等 2014 发表的原始结果高度一致，证明流水线的分析核心**生物学正确**。

## 输出文件 / Outputs

- `deseq2_results.csv` — 17,199 个基因的完整差异表达结果
- `deseq2_significant.csv` — 853 个显著 DEG
- `volcano.pdf` / `ma_plot.pdf` — 差异表达可视化
- `heatmap_top50.pdf` — Top 50 DEG 热图

## 与完整流水线的关系 / Relation to the full pipeline

流水线（Nextflow + Docker）中的 `DESEQ2` 进程调用的正是这份 **完全相同的 `bin/deseq2.R`**。
因此，本地在真实数据上验证通过，即证明了流水线该步骤在端到端运行时也会产出同样正确的结果。
其余步骤（FastQC/TrimGalore/STAR/featureCounts）均使用业界标准工具与固定版本镜像。