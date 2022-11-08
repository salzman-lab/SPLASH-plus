process CONCAT_INTERMEDIARIES {

    // tag "${fastq.simpleName}"
    label "process_low"

    input:
    tuple val(anchor), path(intermediaries)

    output:
    path outfile    , emit: intermediary

    script:
    outfile         = "${anchor}.intermediary"
    """
    rm -rf ${outfile}

    cat *.intermediary > ${outfile}
    """
}
