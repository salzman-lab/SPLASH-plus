process SEQKIT_TRANSLATE {

    label "process_medium"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.2.0--h9ee0642_0':
        'quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0' }"

    input:
    fasta

    output:
    path "*aa"      , emit: aa

    script:
    outfile         = "${fasta}".replaceAll(/_id.fasta/, '.aa')
    """
    seqkit translate \\
        -F \\
        --clean \\
        -f 6 \\
        ${fasta} \\
        > ${outfile}
    """
}
