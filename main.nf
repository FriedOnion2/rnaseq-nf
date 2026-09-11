#!/usr/bin/env nextflow
// =============================================================================
//  rnaseq-nf — 可复现的 RNA-seq 高通量分析流水线
//  FASTQ → QC → 修剪 → 比对(STAR) → 定量(featureCounts) → 差异表达(DESeq2) → 可视化
//
//  Reproducible RNA-seq pipeline:
//  FASTQ → QC → Trim → Align (STAR) → Quantify (featureCounts) → DE (DESeq2) → Viz
// =============================================================================

nextflow.enable.dsl = 2

// -------------------------------------------------------------------------
// 参数定义 / Parameters
// -------------------------------------------------------------------------
params.reads          = "${projectDir}/data/samplesheet.csv"   // 样本表 / samplesheet
params.genome_fasta   = null   // 参考基因组 FASTA / reference genome
params.gtf            = null   // 基因注释 GTF / annotation
params.outdir         = "${launchDir}/results"
params.contrast       = null   // 差异对比，格式 "conditionA,conditionB" 例如 "treat,control"
params.stranded       = 'no'   // 链特异性 / strandedness: no|yes|reverse
params.single_end     = false  // 单端测序 / single-end mode
params.adapter_file   = "${projectDir}/assets/adapters.fa"
params.min_trim_len   = 20
params.star_threads   = 8
params.refdir         = "${launchDir}/refs"   // 存放下载好的参考文件
params.skip_rseqc      = false

include { FASTQC }                    from './modules/fastqc/main'
include { MULTIQC }                   from './modules/multiqc/main'
include { TRIMGALORE }                from './modules/trimgalore/main'
include { STAR_INDEX; STAR_ALIGN }    from './modules/star/main'
include { FEATURECOUNTS }             from './modules/featurecounts/main'
include { DESEQ2 }                    from './modules/deseq2/main'
include { QUALIMAP }                  from './modules/qualimap/main'


workflow {

    // ---- 1. 读入样本表并解析通道 / Read samplesheet, build channels ----
    def raw_samples = Channel
        .fromPath(params.reads, checkIfExists: true)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            def id   = row.sample
            def cond = row.condition
            def fastq1 = row.fastq_1
            def fastq2 = row.fastq_2 ?: null
            [id, cond, fastq1, fastq2]
        }

    def samples = raw_samples
        .map { id, cond, fq1, fq2 -> [id, cond, file(fq1, checkIfExists: true), fq2 ? file(fq2, checkIfExists: true) : null] }

    // ---- 2. 原始数据质控 / Raw QC ----
    FASTQC(samples.map { id, cond, fq1, fq2 -> [id, fq1, fq2] })

    // ---- 3. 修剪接头与低质量碱基 / Trim adapters & low-quality bases ----
    TRIMGALORE(
        samples.map { id, cond, fq1, fq2 -> [id, fq1, fq2] },
        params.adapter_file,
        params.min_trim_len,
        params.single_end
    )

    // ---- 4. 构建 STAR 索引（一次性）/ Build STAR index (once) ----
    // 参考文件：先从 refdir 自动探测，其次用参数显式指定
    def genome_fa = params.genome_fasta
    def annot_gtf = params.gtf
    if (!genome_fa) {
        genome_fa = file("${params.refdir}/genome.fa", checkIfExists: true)
    }
    if (!annot_gtf) {
        annot_gtf = file("${params.refdir}/genes.gtf", checkIfExists: true)
    }
    def ref_ch = Channel.value([genome_fa, annot_gtf])

    STAR_INDEX(ref_ch)

    // ---- 5. 比对 / Align ----
    // 将每个样本与 STAR 索引配对（索引为单值，与 reads 通道 combine 后广播到每个样本）
    STAR_ALIGN(
        TRIMGALORE.out.reads.combine(STAR_INDEX.out.index),
        params.star_threads
    )

    // ---- 6. 定量 / Quantify ----
    // 每个样本与 GTF 配对（GTF 来自 value 通道，自动广播）
    def gtf_ch = ref_ch.map { fa, gtf -> gtf }
    FEATURECOUNTS(
        STAR_ALIGN.out.bam.combine(gtf_ch),
        params.stranded,
        !params.single_end
    )

    // ---- 7. 差异表达分析 / Differential expression ----
    // 收集全部 counts 文件 + 对比参数 + R 脚本
    DESEQ2(
        FEATURECOUNTS.out.counts.collect(),
        params.contrast,
        file("${projectDir}/bin/deseq2.R", checkIfExists: true)
    )

    // ---- 8. 比对后质控 / Post-alignment QC ----
    QUALIMAP(
        STAR_ALIGN.out.bam.combine(gtf_ch)
    )

    // ---- 9. 汇总质控报告 / Aggregate QC report ----
    // 各上游输出为 tuple(sample_id, path)，这里提取纯 path 再汇总
    MULTIQC(
        FASTQC.out.zip.map          { id, html, zip -> zip }.collect(),
        TRIMGALORE.out.reports.map  { id, report -> report }.collect(),
        STAR_ALIGN.out.logs.map     { id, log -> log }.collect(),
        QUALIMAP.out.reports.map    { id, report -> report }.collect()
    )

    // ---- 发布结果 / Publish ----
    MULTIQC.out.report | view
}