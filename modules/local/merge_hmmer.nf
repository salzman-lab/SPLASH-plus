process MERGE_HMMER {

    label "process_low"

    input:
    compactors_fasta
    spindles_fasta
    compactor_summary
    tblouts

    output:
    path outfile        , emit: tsv

    script:
    compactors_tblout   = "compactors.tblout"
    spindles_tblout     = "spindles.tblout"
    outfile             = "compactor_summary_hmmer_merged.tsv"
    """
    pfam_merge.py \\
        --compactor_summary ${compactor_summary} \\
        --compactors_fasta ${compactors_fasta} \\
        --spindles_fasta ${spindles_fasta} \\
        --compactors_tblout ${compactors_tblout} \\
        --spindles_tblout ${spindles_tblout} \\
        --outfile ${outfile}
    """
}
