#!/usr/bin/env bash
# =============================================================================
#  下载演示用真实公开数据 / Download real public demo dataset
#  数据来源：NCBI SRA — 人类细胞系 RNA-seq（treat vs control 小样本子集）
#  使用 SRAtoolkit 的 fasterq-dump 下载并压缩为 gz
#
#  Usage: bash scripts/download_demo_data.sh
# =============================================================================
set -euo pipefail

# 示例：使用一个真实的公开 RNA-seq 数据集（此处以 6 个 SRA 登录号占位）
# 请根据实际要复现的研究替换为对应 SRA accession
SRR_LIST=(
  SRR24231729
  SRR24231730
  SRR24231731
  SRR24231732
  SRR24231733
  SRR24231734
)

mkdir -p data/reads
OUTDIR=data/reads

for srr in "${SRR_LIST[@]}"; do
  echo "下载 $srr ..."
  fasterq-dump "$srr" -O "$OUTDIR" --split-files --threads 4
  gzip "$OUTDIR/${srr}_1.fastq"
  gzip "$OUTDIR/${srr}_2.fastq"
done

echo "完成。使用 data/samplesheet.csv 中的路径映射。"
echo "提示：请将 samplesheet.csv 中的 fastq 路径改为此处生成的文件名。"