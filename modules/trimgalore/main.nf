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
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: reads
    tuple val(sample_id), path("*_trimming_report.txt"), emit: reports
    path "*_trimming_report.txt", emit: report_raw

    script:
    if (single_end) {
        """
        trim_galore \\
            --cores ${task.cpus} \\
            --length ${min_len} \\
            --gzip \\
            --adapter fasta:${adapter_file} \\
            ${reads1}
        """
    } else {
        """
        trim_galore \\
            --cores ${task.cpus} \\
            --length ${min_len} \\
            --gzip \\
            --paired \\
            --retain_unpaired \\
            --adapter AGATCGGAAGAGC \\
            --adapter2 AGATCGGAAGAGC \\
            ${reads1} ${reads2}
        """
    }

    stub:
    """
    touch ${sample_id}_val_1.fq.gz
    touch ${sample_id}_val_2.fq.gz
    touch ${sample_id}_R1_trimming_report.txt
    """
}