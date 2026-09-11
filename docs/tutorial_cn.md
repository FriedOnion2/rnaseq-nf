# 手把手教学：从零复现一个 RNA-seq 差异表达分析流水线

> 本教程假设你**完全没接触过**高通量测序分析。我们从"为什么要这样做"开始，
> 一步步带你跑通、读懂、并改造成你自己的项目。
>
> 目标：搞懂一条 FASTQ → 差异表达基因（DEG）流水线的**每一个环节**。

---

## 目录

- [第 0 课：先建立整体地图](#第-0-课先建立整体地图)
- [第 1 课：生物学问题是"什么变了"](#第-1-课生物学问题是什么变了)
- [第 2 课：环境准备（WSL2 + Docker + Nextflow）](#第-2-课环境准备wsl2--docker--nextflow)
- [第 3 课：准备数据（参考基因组 + 真实 reads）](#第-3-课准备数据参考基因组--真实-reads)
- [第 4 课：一条命令跑通全流程](#第-4-课一条命令跑通全流程)
- [第 5 课：读懂每一步的输出](#第-5-课读懂每一步的输出)
- [第 6 课：读懂代码——main.nf 与模块](#第-6-课读懂代码mainnf-与模块)
- [第 7 课：统计原理——DESeq2 为什么这样做](#第-7-课统计原理deseq2-为什么这样做)
- [第 8 课：换成你自己的数据](#第-8-课换成你自己的数据)
- [第 9 课：复现性与常见坑](#第-9-课复现性与常见坑)

---

## 第 0 课：先建立整体地图

整条流水线回答一个生物学问题：

> **"处理组 vs 对照组，哪些基因的表达量发生了显著变化？"**

数据流一条线：

```
FASTQ(原始序列) → FastQC(质控) → TrimGalore(去接头) → STAR(比对到基因组)
                    → featureCounts(数每个基因几条read) → DESeq2(统计找显著) → 画图
```

先记住 3 个贯穿全程的概念：

| 概念 | 是什么 | 类比 |
|---|---|---|
| **FASTQ** | 测序仪吐出的原始文件，每条 read 4 行（名/序列/+/质量） | 一箱没分类的纸片 |
| **BAM** | 比对结果，记录每条 read 落在基因组哪个位置 | 把纸片都贴到了地图对应位置 |
| **Count 矩阵** | 表：行=基因，列=样本，值=该基因被读到的条数 | 统计每个"邮区"收到几张纸片 |

> **为什么不能跳过"比对"直接数基因？** 因为一条 read 只是一小段约 100–150 个碱基的
> 短序列，我们不知道它来自哪个基因。必须把它"对齐"（align）到参考基因组上，才能知道
> 它落在哪个基因的区间里。

---

## 第 1 课：生物学问题是什么变了

**背景**：本项目用的真实数据是 **GSE52778**（Himes et al. 2014），研究**地塞米松**
（一种糖皮质激素）对人气道平滑肌细胞转录组的影响。样本分两组：

- **control**（对照）：不处理，3 个生物学重复
- **treated**（处理）：加地塞米松，3 个生物学重复

> **为什么每组要 3 个重复？** 生物学个体差异很大。只有一次测量，我们无法区分
> "这个基因真的变了"和"这个样本刚好抽到偏高"。3 个重复是差异表达分析的最低要求，
> 更多重复（6–8）能显著提高统计功效。

**"有含金量"的关键**：地塞米松的转录组效应有大量文献验证过。所以我们可以用
**已知的响应基因**（如 FKBP5、DUSP1、ZBTB16…）作为**外部金标准**，检验我们
流水线跑出来的结果对不对。结果：**9/11 个经典基因全部被判为显著上调** ✅，
证明分析核心生物学正确（详见 `docs/airway_realdata_validation.md`）。

---

## 第 2 课：环境准备（WSL2 + Docker + Nextflow）

本流水线用 Nextflow 编排 + Docker 容器化，**必须跑在 Linux 环境**。

### 2.1 为什么必须是 Linux

- **Nextflow 无法在原生 Windows 运行**（会报 `Unknown signal: HUP`）。
- **Docker 的 Linux 容器后端需要 WSL2**。

### 2.2 装 WSL2（一次性）

管理员 PowerShell 执行，然后重启：

```powershell
wsl --install --no-distribution
wsl --set-default-version 2
```

### 2.3 在 WSL2 里装 Docker

> ⚠️ 关键经验：用 **WSL2 里的原生 Docker Engine**，而不是依赖 Docker Desktop 的
> WSL 集成（后者在无 GUI 时很难配置成功）。

```bash
# 进入 Ubuntu 后
sudo apt-get update
sudo apt-get install -y docker.io
# 启动 Docker 守护进程
sudo dockerd &
# 验证
docker run --rm hello-world   # 应输出 "Hello from Docker!"
```

### 2.4 装 Java + 拿到 Nextflow

```bash
sudo apt-get install -y openjdk-21-jre-headless
chmod +x bin/nextflow   # 仓库已附官方 launcher，首次运行会自举完整发行版
```

> ⚠️ **Nextflow 版本坑**：本项目用 **24.10.8**。Nextflow 26.x 起把 docker 执行器
> 拆成了外部插件，`executor = 'docker'` 会报 `Unknown executor name`。别用 26.x。

---

## 第 3 课：准备数据（参考基因组 + 真实 reads）

流水线需要两类输入：**参考基因组**（拿来比对的"标准地图"）和 **reads**（被测样本）。

### 3.1 下载参考基因组（chr21 子集，几 MB）

```bash
bash scripts/download_reference.sh --mini
```

会得到：
- `refs/genome.fa` — 人类 21 号染色体序列（真实）
- `refs/genes.gtf` — 基因注释，记录每个基因在基因组上的起止位置

> `--mini` 只取 chr21（完整 GRCh38 有 ~3GB，太大）。跑通流程用 chr21 足矣；正式分析
> 去掉 `--mini` 下完整版。

### 3.2 下载真实 reads（子采样 50 万对/样本）

```bash
bash scripts/download_demo_data.sh 500000
```

从 ENA 直连下载 GSE52778 的 6 个样本 FASTQ，每个子采样 50 万条 read（控制体积）。

### 3.3 样本表 `data/samplesheet.csv`

这是流水线的"入口"，告诉它有哪些样本、什么分组、read 在哪：

```csv
sample,condition,fastq_1,fastq_2
control_rep1,control,data/reads/control_rep1_R1.fastq.gz,data/reads/control_rep1_R2.fastq.gz
...
treated_rep1,treated,data/reads/treated_rep1_R1.fastq.gz,data/reads/treated_rep1_R2.fastq.gz
```

> 注意是**双端测序**，所以每个样本有 R1/R2 两个文件（read1 和 read2 是同一段 DNA
> 两端分别测到的）。

---

## 第 4 课：一条命令跑通全流程

```bash
bash scripts/run_pipeline.sh
```

这个脚本做了三件事：
1. 确保 Docker 守护进程在跑
2. 把 work 目录放到 `/nf-work`（Linux 原生分区）
3. 用 `-profile docker -resume` 起 Nextflow

> ⚠️ **两个最关键的坑**（新手必看）：
> 1. **`-work-dir /nf-work`**：STAR 需要创建 FIFO 命名管道，Windows NTFS（`/mnt/d`）
>    不支持 FIFO。所以 work 目录必须放 Linux 原生 ext4 分区（`/nf-work`）。
> 2. **不能用 `/root`**：MultiQC 容器以非 root 用户运行，进不去 700 权限的 `/root`。

手动等价命令（拆开看更清楚）：

```bash
nextflow run main.nf \
  --contrast treated,control \
  -profile docker \
  -work-dir /nf-work \
  -resume
```

### `-resume` 的妙用

Nextflow 会缓存每个任务的输入/命令哈希。你改了什么、它只重跑受影响的部分，其余
从缓存直接取（`cached`）。**改代码后反复调试时，`-resume` 能省几十分钟。**

---

## 第 5 课：读懂每一步的输出

运行完，`results/` 长这样（每步一个子目录）：

```
results/
├── deseq2/          # 核心产出
│   ├── deseq2_results.csv      # 所有基因的差异表达表
│   ├── deseq2_significant.csv  # 显著 DEG
│   ├── volcano.pdf             # 火山图
│   ├── ma_plot.pdf             # MA 图
│   └── heatmap_top50.pdf       # Top50 热图
├── fastqc/          # 原始数据质量
├── trimgalore/      # 修剪后数据 + 报告
├── star/            # BAM + 比对率日志（Log.final.out）
├── featurecounts/   # 基因 counts 矩阵
├── qualimap/        # 比对后质控
└── multiqc/multiqc_report.html  # 全流程 QC 汇总（先看这个！）
```

### 逐个读懂关键文件

**（1）`deseq2_results.csv` 是最终答案**，每行一个基因，关键列：

| 列 | 含义 | 怎么读 |
|---|---|---|
| `log2FoldChange` | 处理组相对对照组的表达倍数，取 log2 | +1 = 上调 2 倍，-1 = 下调 2 倍 |
| `padj` | 校正后 p 值 | < 0.05 才叫"显著" |
| `baseMean` | 归一化后平均表达量 | 越大表达越高 |

本次演示结果（chr21 子集）检出 4 个显著 DEG：

```
ENSG00000154734  log2FC=+2.31  padj=1.03e-17   ⬆ 上调约 5 倍
ENSG00000197381  log2FC=+2.33  padj=4.69e-04   ⬆ 上调
ENSG00000237569  log2FC=-1.18  padj=6.06e-03   ⬇ 下调
ENSG00000159200  log2FC=+1.30  padj=4.14e-02   ⬆ 上调
```

> 为什么完整数据的验证有 853 个 DEG，而这次演示只有 4 个？因为这次为了快速跑通，
> 只用了 chr21（约全基因组 1/50 的基因）+ 子采样 reads。**链路正确性**已经证明；
> 数量少是"数据缩小了"，不是"方法错了"。

**（2）`star/*.Log.final.out`** — 看比对率。里面 `Uniquely mapped reads %` 通常在
80%+ 才算好数据。

**（3）`multiqc_report.html`** — 用浏览器打开，所有 QC 一页汇总，先看有没有红色警告。

---

## 第 6 课：读懂代码——main.nf 与模块

Nextflow 用 `main.nf` 定义"流程编排"，每个工具一个 `modules/<工具>/main.nf` 模块。

### 6.1 `main.nf` — 流水线的"骨架"

```groovy
workflow {
    // 1. 读样本表 → 通道
    def samples = Channel.fromPath(params.reads).splitCsv(...).map {...}

    // 2-9. 把通道接起来（每个模块是一个"加工站"）
    FASTQC(samples...)
    TRIMGALORE(samples...)
    STAR_INDEX(ref_ch)
    STAR_ALIGN(TRIMGALORE.out.reads.combine(STAR_INDEX.out.index), ...)
    FEATURECOUNTS(STAR_ALIGN.out.bam.combine(gtf_ch), ...)
    DESEQ2(FEATURECOUNTS.out.counts.collect(), params.contrast, ...)
    QUALIMAP(STAR_ALIGN.out.bam.combine(gtf_ch))
    MULTIQC(...)
}
```

**核心思想（数据流）**：每个模块的**输出**就是下一个模块的**输入**。Nextflow 用
**通道（Channel）** 传递数据，且能**自动并行**（6 个样本的 FastQC 会同时跑）。

> 关键点 `combine`：STAR 索引是"一个"，但有"6 个样本"，`combine` 把单值索引
> 广播到每个样本——"一个索引 + 每个 reads = 配对"。

### 6.2 一个模块长什么样（以 STAR_ALIGN 为例）

```groovy
process STAR_ALIGN {
    tag "$sample_id"                       // 日志里显示样本名
    container '.../star:2.7.11b...'        // 用哪个 Docker 镜像
    publishDir "${params.outdir}/star"     // 结果发布到哪里

    input:  tuple val(sample_id), path(reads1), path(reads2), path(star_index)
    output: tuple val(sample_id), path("*.Aligned.sortedByCoord.out.bam")

    script: """
    STAR --genomeDir ${star_index} --readFilesIn ${reads1} ${reads2} ...
    """
}
```

三个核心块：
- **`input`**：声明吃什么（类型：`val`=值、`path`=文件）
- **`output`**：声明吐什么（glob 通配符匹配实际生成的文件名）
- **`script`**：真正执行的 shell 命令（在容器里跑）

> 新手最常踩的坑就是 **output 的 glob 和实际文件名不匹配**（比如 `*_Log.final.out`
> vs 实际 `*.Log.final.out`）——文件名对不上会报 `Missing output file`。

### 6.3 为什么用 Docker 容器

每个模块固定一个镜像 + 版本号（如 `star:2.7.11b--h5ca1c30_4`），保证
**任何人、任何机器、任何时间**跑出来的结果一致。这就是"可复现"的核心。

### 6.4 完整逐行注释 `main.nf`

下面把本仓库真实的 `main.nf` 完整拆开，一句一句讲。读完这一节，你就真正"看懂"了流水线。

```groovy
nextflow.enable.dsl = 2
// ↑ 声明用 DSL2 语法（现代 Nextflow 标准）

params.reads    = "${projectDir}/data/samplesheet.csv"
params.contrast = null          // 例如 "treated,control"
params.refdir   = "${launchDir}/refs"
// ↑ params = 可被命令行 --xxx 覆盖的全局参数

include { FASTQC }      from './modules/fastqc/main'
include { STAR_INDEX; STAR_ALIGN } from './modules/star/main'
// ↑ include = 把各模块"注册"进来。STAR 模块里同时有 STAR_INDEX 和 STAR_ALIGN 两个 process

workflow {
    // ---- ① 读样本表 → 通道 ----
    def raw_samples = Channel
        .fromPath(params.reads, checkIfExists: true)   // 找到 CSV 文件
        .splitCsv(header: true, sep: ',', strip: true) // 按行解析成记录
        .map { row -> [row.sample, row.condition, row.fastq_1, row.fastq_2] }
    // ↑ 每一行 CSV 变成一个 [样本名, 条件, R1, R2] 四元组，流进通道

    def samples = raw_samples.map { id, cond, fq1, fq2 ->
        [id, cond, file(fq1), fq2 ? file(fq2) : null]
    }
    // ↑ file() 把字符串路径变成"文件对象"，Nextflow 才会自动追踪/传递真实文件

    // ---- ② 原始质控 ----
    FASTQC(samples.map { id, cond, fq1, fq2 -> [id, fq1, fq2] })
    // ↑ map 把 [id, cond, fq1, fq2] 改成 [id, fq1, fq2]（FastQC 不需要 cond）

    // ---- ③ 修剪 ----
    TRIMGALORE(
        samples.map { id, cond, fq1, fq2 -> [id, fq1, fq2] },
        params.adapter_file, params.min_trim_len, params.single_end
    )

    // ---- ④ 建 STAR 索引（只跑一次）----
    def genome_fa = params.genome_fasta ?: file("${params.refdir}/genome.fa")
    def annot_gtf = params.gtf         ?: file("${params.refdir}/genes.gtf")
    def ref_ch = Channel.value([genome_fa, annot_gtf])
    // ↑ Channel.value = "值通道"：内容固定只有一个，供所有样本复用
    STAR_INDEX(ref_ch)

    // ---- ⑤ 比对 ----
    STAR_ALIGN(
        TRIMGALORE.out.reads.combine(STAR_INDEX.out.index),
        params.star_threads
    )
    // ↑ combine：把"6 个样本的 reads" × "1 个索引" → 6 个配对。
    //   这是关键技巧：单值索引被"广播"到每个样本。

    // ---- ⑥ 定量 ----
    def gtf_ch = ref_ch.map { fa, gtf -> gtf }   // 从值通道里只取 GTF
    FEATURECOUNTS(
        STAR_ALIGN.out.bam.combine(gtf_ch),
        params.stranded, !params.single_end
    )

    // ---- ⑦ 差异表达 ----
    DESEQ2(
        FEATURECOUNTS.out.counts.collect(),   // collect：把 6 个 counts 聚成一个列表
        params.contrast,                       // "treated,control"
        file("${projectDir}/bin/deseq2.R")     // R 脚本
    )
    // ↑ 注意 DESeq2 和前面不同：它要"所有样本一起"分析，所以用 collect()

    // ---- ⑧ 比对后质控 ----
    QUALIMAP(STAR_ALIGN.out.bam.combine(gtf_ch))

    // ---- ⑨ 汇总质控 ----
    MULTIQC(
        FASTQC.out.zip.map       { id, html, zip -> zip }.collect(),
        TRIMGALORE.out.reports.map { id, r -> r }.collect(),
        STAR_ALIGN.out.logs.map   { id, log -> log }.collect(),
        QUALIMAP.out.reports.map  { id, r -> r }.collect()
    )
    // ↑ 每个上游输出都是 tuple(id, path)，这里 map 只留 path，再 collect 汇总
}
```

**三个贯穿全程的 Nextflow 概念：**

| 概念 | 作用 | 本流水线例子 |
|---|---|---|
| **Channel（通道）** | 数据在模块间流动的管道 | `samples`、`ref_ch` |
| **`.map{}`** | 逐个变换通道里的元素 | `[id,cond,fq1,fq2]` → `[id,fq1,fq2]` |
| **`.combine()`** | 两个通道做笛卡尔配对 | 6 样本 × 1 索引 |
| **`.collect()`** | 把流里所有元素攒成一个列表 | 6 个 counts → 1 个列表给 DESeq2 |

> **为什么 `combine` 和 `collect` 是两个不同动作？**
> - `combine`：每个样本要和索引/GTF **一对一**配对（保持 6 份独立）
> - `collect`：DESeq2 需要**一次性拿到全部 6 个 counts**（合并分析）
>
> 理解这两个的区别，就理解了 Nextflow 数据流的精髓。

---

## 第 7 课：统计原理——DESeq2 为什么这样做

这是整条流水线"最烧脑"但也"最有含金量"的一步。

### 7.1 为什么不能直接比较 count 大小？

假设基因 A 在对照样本里数到 1000 条，处理样本里数到 2000 条——能说它上调 2 倍吗？

**不能。** 因为两个样本**测序深度**可能不同（一个测了 1 千万条，一个测了 2 千万条）。
所以要先**归一化**：把每个样本的总量拉到同一水平（DESeq2 用"中位比值法"）。

### 7.2 为什么用"负二项分布"而不是正态分布？

RNA-seq 的 count 数据有两个特征：
1. **离散**（非负整数，不能是 3.7 条）
2. **方差大于均值**（低表达基因波动特别大，"过离散" overdispersion）

所以用**负二项分布（Negative Binomial）**建模，它有两个参数：
- 均值 μ（基因的平均表达）
- 离散度 dispersion（表达越低的基因，dispersion 越大，越不稳定）

> 直觉理解：正态分布假设"波动是均匀的"，但低表达基因 count 可能是 0, 3, 50 这样
> 剧烈跳动，负二项分布能刻画这种"越穷越不稳定"的特性。

### 7.3 p 值为什么要校正（padj）？

我们一次测了 1.7 万个基因，每个都算一个 p 值。即使**完全没有差异**，按 0.05 阈值，
也会"碰巧"有 ~850 个基因的 p < 0.05（假阳性）。

**多重检验校正**（Benjamini-Hochberg，FDR）就是控制这批假阳性的比例。所以看结果
一定要看 **`padj`（校正后 p 值），而不是原始 `pvalue`**。

### 7.4 对比的方向（contrast）

`--contrast treated,control` 的语义是 **"treated 相对 control"**，即：

```
log2FC > 0  → 处理组上调
log2FC < 0  → 处理组下调
```

顺序千万别搞反，否则上调和下调会整体颠倒。

---

## 第 8 课：换成你自己的数据

### 8.1 只改三个地方

1. **准备你自己的 FASTQ**，放到 `data/reads/`
2. **改 `data/samplesheet.csv`**，填你的样本名、分组、文件路径
3. **改 `--contrast`** 参数，换成你的两个组名（例如 `--contrast mutant,wildtype`）

### 8.2 换参考基因组（如果你不是人）

```bash
# 下载你自己的物种 FASTA + GTF，放到 refs/genome.fa 和 refs/genes.gtf
# 或用完整 GRCh38（去掉 --mini）
bash scripts/download_reference.sh   # 不带 --mini 就是完整人类基因组
```

### 8.3 常见调整

| 想改什么 | 在哪改 |
|---|---|
| 单端测序 | `--single_end` 传 true |
| 链特异性文库 | `--stranded reverse` |
| 显著阈值（默认 padj<0.05） | `bin/deseq2.R` 里的 `0.05` |
| 比对线程数 | `--star_threads` |

---

## 第 9 课：复现性与常见坑

### 复现性三要素（本项目都做到了）

1. **固定版本**：每个工具镜像都钉死版本号
2. **容器化**：环境差异被 Docker 隔离
3. **记录参数**：`results/pipeline_info/execution_report.html` 记录每步命令

### 高频坑速查

| 报错 | 原因 | 解决 |
|---|---|---|
| `could not create FIFO file` | work 目录在 NTFS | `-work-dir /nf-work` |
| `Unknown executor name: docker` | 用了 Nextflow 26.x | 换回 24.10.8 |
| FastQC 卡死 `No fonts found` | biocontainers/fastqc 缺字体 | 用 `staphb/fastqc` |
| `.command.run: Permission denied` | work 在 `/root` | 放 `/nf-work` 并 chmod 777 |
| 镜像 `not found` | tag 失效 | `docker manifest inspect` 查新 tag |
| `Missing output file` | output glob 和实际文件名不匹配 | 核对文件名 |

---

## 下一步 / 深入阅读

- [`docs/airway_realdata_validation.md`](airway_realdata_validation.md) — 真实数据的生物学金标准验证
- [`docs/local_deseq2_verification.md`](local_deseq2_verification.md) — 合成数据的零假阳性验证
- [DESeq2 原始论文](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) — Love et al. 2014
- [Nextflow 官方文档](https://www.nextflow.io/docs/latest/) — DSL2 与通道