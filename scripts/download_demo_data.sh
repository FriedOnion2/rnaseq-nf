#!/usr/bin/env bash
# =============================================================================
#  下载真实公开 RNA-seq 演示数据（GSE52778）—— 直接用 ENA API 查询真实路径
#  API 返回两列（run_accession / fastq_ftp，tab 分隔）；fastq_ftp 内多个文件用 ; 分隔，
#  我们取双端配对文件 _1 与 _2。
#
#  数据来源：Himes et al. 2014, PLoS One — GEO GSE52778
#  Usage: bash scripts/download_demo_data.sh [每样本read对数]
# =============================================================================
set -uo pipefail

PAIRS="${1:-500000}"
LINES=$((PAIRS * 4))
OUTDIR="data/reads"
mkdir -p "$OUTDIR"

declare -A RUN2SAMPLE=(
  [SRR1039508]=control_rep1
  [SRR1039509]=treated_rep1
  [SRR1039512]=control_rep2
  [SRR1039513]=treated_rep2
  [SRR1039516]=control_rep3
  [SRR1039517]=treated_rep3
)

for acc in "${!RUN2SAMPLE[@]}"; do
  sample="${RUN2SAMPLE[$acc]}"
  echo "[下载] ${acc} -> ${sample}"

  # 取 fastq_ftp 列（第 2 列），去掉可能的引号/回车
  ftp=$(curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp" \
        | awk -F'\t' 'NR==2{print $2}' | tr -d '\r')

  # 从 ; 分隔的列表里选 _1 和 _2 文件
  r1=""; r2=""
  IFS=';' read -ra parts <<< "$ftp"
  for p in "${parts[@]}"; do
    case "$p" in
      *_1.fastq.gz) r1="https://$p" ;;
      *_2.fastq.gz) r2="https://$p" ;;
    esac
  done

  if [[ -z "$r1" || -z "$r2" ]]; then
    echo "  !! 未找到双端配对 _1/_2 文件，跳过。ftplist=$ftp"
    continue
  fi

  echo "  R1: $r1"
  echo "  R2: $r2"

  curl -fsSL "$r1" | zcat | head -n "$LINES" | gzip > "${OUTDIR}/${sample}_R1.fastq.gz"
  curl -fsSL "$r2" | zcat | head -n "$LINES" | gzip > "${OUTDIR}/${sample}_R2.fastq.gz"
  echo "  -> 完成 ($(du -h "${OUTDIR}/${sample}_R1.fastq.gz" | cut -f1) / R2)"
done

echo ""
echo "全部下载完成："
ls -la "$OUTDIR"