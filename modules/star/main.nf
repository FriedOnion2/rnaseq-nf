// =============================================================================
//  STAR 模块 — 索引构建 + 比对 / Genome index build + alignment
// =============================================================================

// ---- STAR 索引：一次性构建，供所有样本复用 ----
process STAR_INDEX {
    label 'process_high'
    label 'process_long'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_4'

    input:
    tuple path(genome_fasta), path(gtf)

    output:
    path "star_index/", emit: index

    script:
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus} \\
        --sjdbOverhang 99 \\
        --genomeSAindexNbases 12
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/SA star_index/Genome
    """
}

// ---- STAR 比对 / Alignment ----
process STAR_ALIGN {
    tag "$sample_id"
    label 'process_high'
    label 'process_long'

    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_4'

    input:
    tuple val(sample_id), path(reads1), path(reads2), path(star_index)
    val   threads

    output:
    tuple val(sample_id), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample_id), path("*_Log.final.out"), emit: logs
    path "*_Log.final.out", emit: log_raw
    path "*.Aligned.sortedByCoord.out.bam", emit: bam_raw

    script:
    def fq = reads2 ? "${reads1} ${reads2}" : "${reads1}"
    def readtype = reads2 ? '--readFilesIn' : '--readFilesIn'
    """
    STAR \\
        --genomeDir ${star_index} \\
        ${readtype} ${fq} \\
        --readFilesCommand zcat \\
        --runThreadN ${threads} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD \\
        --quantMode TranscriptomeSAM GeneCounts \\
        --outFileNamePrefix ${sample_id}. \\
        --outReadsUnmapped Fastx
    """

    stub:
    """
    touch ${sample_id}.Aligned.sortedByCoord.out.bam
    touch ${sample_id}_Log.final.out
    """
}