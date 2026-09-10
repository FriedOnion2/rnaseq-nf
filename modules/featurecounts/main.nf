// =============================================================================
//  featureCounts 模块 — 基因级别定量 / Gene-level read quantification
// =============================================================================

process FEATURECOUNTS {
    tag "$sample_id"
    label 'process_medium'

    container 'quay.io/subread_sourceforge/subread:2.0.6'

    input:
    tuple val(sample_id), path(bam), path(gtf)
    val   stranded
    val   paired

    output:
    tuple val(sample_id), path("*_counts.txt"), emit: counts_raw
    path "*counts.txt", emit: counts
    path "*counts.txt.summary", emit: summary

    script:
    def strand_flag = '0'
    if (stranded == 'yes') { strand_flag = '1' }
    else if (stranded == 'reverse') { strand_flag = '2' }
    def pair_flag = paired ? '-p --countReadPairs' : ''
    """
    featureCounts \\
        -a ${gtf} \\
        -o ${sample_id}_counts.txt \\
        -t exon \\
        -g gene_id \\
        -s ${strand_flag} \\
        ${pair_flag} \\
        -T ${task.cpus} \\
        ${bam}
    """

    stub:
    """
    touch ${sample_id}_counts.txt
    touch ${sample_id}_counts.txt.summary
    """
}