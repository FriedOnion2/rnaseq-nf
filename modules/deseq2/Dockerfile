# Dockerfile — 为 DESeq2 步骤提供自定义镜像
# 将仓库内的 deseq2.R 打进镜像，使模块可直接调用 Rscript /opt/deseq2.R
#
# Build:  docker build -f modules/deseq2/Dockerfile -t rnaseq-nf/deseq2:1.0 .
FROM quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hdfd78af_0

LABEL org.opencontainers.image.title="rnaseq-nf DESeq2"
LABEL org.opencontainers.image.description="DESeq2 differential expression for rnaseq-nf"

# 安装 R 依赖（optparse 等）—— 基础镜像可能未包含
RUN Rscript -e 'install.packages("optparse", repos="https://cloud.r-project.org")' \
    || echo "optparse 可能已存在，跳过"

# 拷贝分析脚本
COPY bin/deseq2.R /opt/deseq2.R

WORKDIR /workspace

ENTRYPOINT ["Rscript", "/opt/deseq2.R"]