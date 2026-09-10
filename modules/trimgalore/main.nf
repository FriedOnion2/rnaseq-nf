// =============================================================================
//  TrimGalore 模块 — 修剪接头与低质量碱基 / Trim adapters + low-quality bases
// =============================================================================

process TRIMGALORE {
    tag "$sample_id"
    label 'process_medium'

    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_4'

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path  adapter_file
    val   min_len
    val   single_end

    output:
    tuple val(sample_id), path("*val*.fq.gz"), path("*_2_val*.fq.gz"), emit: reads
    tuple val(sample_id), path("*_trimming_report.txt"), emit: reports
    path "*_trimming_report.txt", emit: report_raw
    path "*.fq.gz", emit: fastq

    script:
    def pe = single_end ? '' : '--paired'
    def fq2 = (!single_end && reads2) ? "${reads2}" : ''
    if (single_end) {
        """
        trim_galore \\
            --cores ${task.cpus} \\
            --length ${min_len} \\
            --gzip \\
            --fastqc \\
            -a file:${adapter_file} \\
            ${reads1}
        """
    } else {
        """
        trim_galore \\
            --cores ${task.cpus} \\
            --length ${min_len} \\
            --gzip \\
            --fastqc \\
            --paired \\
            -a file:${adapter_file} \\
            -a2 file:${adapter_file} \\
            ${reads1} ${fq2}
        """
    }
}