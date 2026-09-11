#!/usr/bin/env bash
# =============================================================================
#  一键启动本流水线（在 WSL2/Ubuntu 内执行）
#  自动确保 dockerd 运行，然后跑 Nextflow 全流程。
#
#  Usage: bash scripts/run_pipeline.sh
# =============================================================================
set -uo pipefail
cd /mnt/d/fuxian/rnaseq-nf || { echo "无法进入项目目录"; exit 1; }

# 1) 确保 dockerd 运行
if ! docker info >/dev/null 2>&1; then
  echo "=== 启动 dockerd ==="
  nohup dockerd > /var/log/dockerd.log 2>&1 &
  for i in $(seq 1 30); do
    docker info >/dev/null 2>&1 && break
    sleep 1
  done
fi
echo "=== Docker 就绪: $(docker info --format '{{.ServerVersion}}') ==="

# 2) Nextflow（bin/nextflow 是官方 launcher 脚本，首次运行会自举完整发行版）
NF_BIN="bin/nextflow"
chmod +x "$NF_BIN"

# 关键：work 目录必须放在 Linux 原生 ext4 分区（而非 /mnt/d 的 NTFS），
# 因为 STAR 需要创建 FIFO 命名管道，NTFS 不支持 FIFO。
# 且必须放在所有用户可访问的路径（不能是 /root，其为 700 权限）：
# 部分镜像（如 MultiQC）以非 root 用户运行，无法进入 /root。
NATIVE_WORK="/nf-work"
mkdir -p "$NATIVE_WORK"
chmod 777 "$NATIVE_WORK"

# 3) 运行流水线（参考文件自动从 refs/ 探测；对比 treated vs control）
echo "=== 启动 Nextflow 流水线 ==="
./"$NF_BIN" run main.nf \
  --contrast treated,control \
  -profile docker \
  -work-dir "$NATIVE_WORK" \
  -resume

echo "=== 完成。结果在 results/ ==="