# =============================================================================
#  环境准备脚本（Windows 专用）/ One-time environment setup (Windows)
#  在 PowerShell（管理员）中运行。
#
#  作用：安装 WSL2 后端。之后：
#    1. 重启电脑
#    2. 启动 Docker Desktop（设置里启用 WSL2 后端）
#    3. 所有 Nextflow 命令都在 WSL2 终端里执行（Nextflow 无法在原生 Windows 运行）
# =============================================================================

Write-Host "=== 1/3 安装 WSL2（需要重启）==="
wsl --install --no-distribution

Write-Host ""
Write-Host "=== 2/3 设置 WSL2 为默认版本 ==="
wsl --set-default-version 2

Write-Host ""
Write-Host "=== 3/3 提示 ==="
Write-Host "请按顺序操作："
Write-Host "  1) 重启电脑"
Write-Host "  2) 启动 Docker Desktop，设置 -> General -> 使用 WSL 2 based engine"
Write-Host "  3) 打开 WSL2 终端 (wsl)，验证:"
Write-Host "       docker run hello-world"
Write-Host "     应输出 'Hello from Docker!'"
Write-Host ""
Write-Host "  4) 在 WSL2 里安装 Nextflow 并跑流水线:"
Write-Host "       cd /mnt/d/fuxian/rnaseq-nf"
Write-Host "       bash scripts/install_nextflow.sh"
Write-Host "       bash scripts/download_reference.sh --mini"
Write-Host "       ./bin/nextflow run main.nf --genome_fasta refs/genome.fa --gtf refs/genes.gtf --contrast treat,control -profile docker"
Write-Host ""
Write-Host "  （快速校验接线/语法，无需 Docker）:"
Write-Host "       ./bin/nextflow run main.nf -profile test -stub-run"