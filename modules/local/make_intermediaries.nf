process MAKE_INTERMEDIARIES {

    tag "${fastq.simpleName}"
    label "process_medium"

    input:
    tuple val(fastq_id), path(fastq)
    path anchors_file
    val anchor_length
    val num_reads
    val max_anchor_reads

    output:
    path "*intermediary"    , emit: intermediary
    path outfile_log        , emit: log

    script:
    outfile_log             = "make_intermediaries_${fastq_id}.log"
    """
    make_intermediaries.py \\
        --fastq ${fastq} \\
        --fastq_id ${fastq_id} \\
        --anchors_file ${anchors_file} \\
        --anchor_length ${anchor_length} \\
        --num_reads ${num_reads} \\
        --max_anchor_reads ${max_anchor_reads} \\
        --outfile_log ${outfile_log}
    """
}
