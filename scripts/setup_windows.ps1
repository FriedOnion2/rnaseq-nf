#!/usr/bin/env bash
# =============================================================================
#  环境准备脚本（Windows 专用）/ One-time environment setup (Windows)
#  在 PowerShell（管理员）中运行：
#    powershell -Command "Start-Process powershell -Verb RunAs -ArgumentList 'scripts/setup_windows.ps1'"
#  或用管理员权限直接运行本脚本。
#
#  作用：安装 WSL2 后端并启动 Docker Desktop（Nextflow 容器运行的前提）。
# =============================================================================

Write-Host "=== 1/3 安装 WSL2（需要重启）==="
wsl --install --no-distribution

Write-Host ""
Write-Host "=== 2/3 设置 WSL2 为默认版本 ==="
wsl --set-default-version 2

Write-Host ""
Write-Host "=== 3/3 提示 ==="
Write-Host "请重启电脑，然后启动 Docker Desktop，等待托盘图标变绿。"
Write-Host ""
Write-Host "验证（重启后）:"
Write-Host "  docker run hello-world"
Write-Host "  应该输出 'Hello from Docker!'"
Write-Host ""
Write-Host "然后即可运行流水线:"
Write-Host "  cd D:\fuxian\rnaseq-nf"
Write-Host "  java -jar bin\nextflow run main.nf --genome_fasta refs/genome.fa --gtf refs/genes.gtf --contrast treat,control -profile docker"