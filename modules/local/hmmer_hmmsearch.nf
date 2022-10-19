process HMMER_HMMSEARCH {

    label "process_high"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1' :
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    aa
    hmmfile

    output:
    path "${tblout}"    , emit: tblout
    path "${stdout}"    , emit: stdout

    script:
    tblout              = "${aa.simpleName}.tblout"
    stdout              = "${aa.simpleName}.stdout"
    """
    hmmsearch \\
        --notextw \\
        -o ${stdout} \\
        --tblout ${tblout} \\
        --cpu ${task.cpus} \\
        ${hmmfile} \\
        ${aa}
    """
}
