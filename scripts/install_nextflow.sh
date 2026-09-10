#!/usr/bin/env bash
# =============================================================================
#  安装 Nextflow 到 bin/nextflow
#  用法: bash scripts/install_nextflow.sh [版本]
#  说明: Nextflow 需要 Java 11–21。脚本从 GitHub Release 下载官方 all-in-one JAR。
# =============================================================================
set -euo pipefail

NF_VERSION="${1:-26.04.6}"
BIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../bin" && pwd)"
URL="https://github.com/nextflow-io/nextflow/releases/download/v${NF_VERSION}/nextflow-${NF_VERSION}-dist"

echo "下载 Nextflow v${NF_VERSION} ..."
curl -sL --retry 5 --retry-delay 2 -o "${BIN_DIR}/nextflow" "${URL}"

chmod +x "${BIN_DIR}/nextflow"
echo "完成：${BIN_DIR}/nextflow"

# 验证
java -version 2>&1 | head -n 1
echo "运行方式:"
echo "  java -jar ${BIN_DIR}/nextflow run main.nf ..."
echo "  或（若已将 bin/ 加入 PATH）: nextflow run main.nf ..."