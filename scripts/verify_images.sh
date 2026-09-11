#!/usr/bin/env bash
# 验证各模块引用的容器镜像 tag 是否真实存在（docker manifest inspect）
set -u

images=(
  "quay.io/biocontainers/fastqc:0.11.6--pl5.22.0_1"
  "quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_2"
  "quay.io/biocontainers/star:2.7.11b--h5ca1c30_4"
  "quay.io/biocontainers/subread:2.0.6--h577a1d6_3"
  "quay.io/biocontainers/bioconductor-deseq2:1.46.0--r44he5774e6_1"
  "quay.io/biocontainers/qualimap:2.3--hdfd78af_0"
  "multiqc/multiqc:v1.25.1"
)

for img in "${images[@]}"; do
  if docker manifest inspect "$img" >/dev/null 2>&1; then
    echo "OK     $img"
  else
    echo "MISSING $img"
  fi
done