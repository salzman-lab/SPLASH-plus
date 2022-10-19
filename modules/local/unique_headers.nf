process UNIQUE_HEADERS {

    label "process_low"

    input:
    compactors_fasta
    spindles_fasta

    output:
    path "*_id.fasta"    , emit: fasta

    script:
    """
    make_unique_headers.py \\
        --compactors_fasta ${compactors_fasta} \\
        --spindles_fasta ${spindles_fasta}
    """
}
