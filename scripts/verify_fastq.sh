#!/usr/bin/env bash
# 临时验证脚本：检查 FASTQ 完整性与 read 数
set -u
cd /mnt/d/fuxian/rnaseq-nf || exit 1
echo "=== FASTQ 完整性验证 ==="
for f in data/reads/*.fastq.gz; do
  lines=$(zcat "$f" 2>/dev/null | wc -l)
  reads=$((lines / 4))
  printf "%-40s %8d reads\n" "$(basename "$f")" "$reads"
done