#!/usr/bin/env bash
# =============================================================================
#  参考基因组下载脚本 / Reference genome download
#  下载 GRCh38 参考 FASTA + 注释 GTF（用于演示，选择 chr21 小片段以节省时间/空间）
#
#  Usage:
#    bash scripts/download_reference.sh        # 完整 GRCh38（约 3GB+，耗时）
#    bash scripts/download_reference.sh --mini  # 演示用 chr21 子集（推荐，几 MB）
# =============================================================================
set -euo pipefail

ENSEMBL_BASE="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna"
GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

MINI_FA="Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
FULL_FA="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

OUTDIR="refs"
MINI=false
for arg in "$@"; do
  case "$arg" in
    --mini) MINI=true ;;
    -o|--outdir) shift; OUTDIR="${1:-refs}" ;;
  esac
done

mkdir -p "$OUTDIR"

echo "[1/3] 下载参考基因组 FASTA ..."
if $MINI; then
  wget -c -q --show-progress -O "$OUTDIR/chr21.fa.gz" "$ENSEMBL_BASE/$MINI_FA"
  gunzip -f "$OUTDIR/chr21.fa.gz"
  mv "$OUTDIR/chr21.fa" "$OUTDIR/genome.fa"
else
  wget -c -q --show-progress -O "$OUTDIR/genome.fa.gz" "$ENSEMBL_BASE/$FULL_FA"
  gunzip -f "$OUTDIR/genome.fa.gz"
fi

echo "[2/3] 下载基因注释 GTF ..."
wget -c -q --show-progress -O "$OUTDIR/Homo_sapiens.GRCh38.110.gtf.gz" "$GTF_URL"
gunzip -f "$OUTDIR/Homo_sapiens.GRCh38.110.gtf.gz"

# 若为 mini 模式，过滤 GTF 只保留 chr21
if $MINI; then
  echo "[2.5/3] 过滤注释至 chr21 ..."
  grep -P '^21\t' "$OUTDIR/Homo_sapiens.GRCh38.110.gtf" > "$OUTDIR/genes_chr21.gtf"
  mv "$OUTDIR/genes_chr21.gtf" "$OUTDIR/genes.gtf"
else
  mv "$OUTDIR/Homo_sapiens.GRCh38.110.gtf" "$OUTDIR/genes.gtf"
fi

echo "[3/3] 完成。参考文件位于 $OUTDIR/"