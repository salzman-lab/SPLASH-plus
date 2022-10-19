include { MAKE_INTERMEDIARIES } from '../../modules/local/make_intermediaries'
include { GENERATE_COMPACTORS } from '../../modules/local/generate_compactors'

workflow GET_COMPACTORS {
    take:
    ch_fastqs

    main:
    /*
    // Process: For each FASTQ file, make intermediary files
    */
    MAKE_INTERMEDIARIES(
        ch_fastqs,
        params.anchors_file,
        params.anchor_length,
        params.num_reads,
        params.max_anchor_reads
    )

    // For each anchor,
    // merge intermediary files from all FASTQ files into one intermediary file per anchor
    ch_intermediaries = MAKE_INTERMEDIARIES.out.intermediary
        .flatten()
        .map{ intermediary ->
            def anchor = intermediary.simpleName.toString().tokenize('_')[0]
            return tuple(anchor, file(intermediary))
        }
        .collectFile(keepHeader: false, storeDir: "${params.outdir}/intermediary"){
            anchor, intermediary ->
            ["${anchor}.intermediary", intermediary]
        }
        .buffer(
            size: 80,
            remainder: true
        )

    /*
    // Process: For each intermediary file, generate compactors
    */
    GENERATE_COMPACTORS(
        ch_intermediaries,
        file(params.fastq_samplesheet),
        params.anchor_length,
        params.epsilon,
        params.branch_threshold,
        params.recursive_depth
    )

    // Merge all compactor output files
    sample_specificity_file = GENERATE_COMPACTORS.out.sample_specificity
        .collectFile(
            name: "${params.outdir}/sample_specificity.tsv",
            keepHeader: true,
            skip: 1
        )

    compactor_summary_file = GENERATE_COMPACTORS.out.compactor_summary
        .collectFile(
            name: "${params.outdir}/compactor_summary.tsv",
            keepHeader: true,
            skip: 1
        )

    compactors_fasta_file = GENERATE_COMPACTORS.out.compactors_fasta
        .collectFile(
            name: "compactor.fasta",
            storeDir: "${params.outdir}"
        )

    spindles_fasta_file = GENERATE_COMPACTORS.out.spindles_fasta
        .collectFile(
            name: "spindles.fasta",
            storeDir: "${params.outdir}"
        )

    emit:
    sample_specificity  = sample_specificity_file
    compactor_summary   = compactor_summary_file
    compactors_fasta    = compactors_fasta_file
    spindles_fasta      = spindles_fasta_file

}

