process GENERATE_COMPACTORS {

    label "process_medium"

    input:
    path intermediaries
    path fastq_samplesheet
    val anchor_length
    val epsilon
    val branch_threshold
    val recursive_depth

    output:
    path "sample_specificity.tsv"   , emit: sample_specificity  , optional: true
    path "compactor.fasta"          , emit: compactors_fasta    , optional: true
    path "spindles.fasta"           , emit: spindles_fasta      , optional: true
    path "compactor_summary.tsv"    , emit: compactor_summary   , optional: true
    path "*.compactor"              , emit: compactor           , optional: true

    script:
    """
    generate_compactors.py \\
        --fastq_samplesheet ${fastq_samplesheet} \\
        --anchor_length ${anchor_length} \\
        --epsilon ${epsilon} \\
        --branch_threshold ${branch_threshold} \\
        --recursive_depth ${recursive_depth}
    """
}
