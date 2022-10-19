include { UNIQUE_HEADERS    } from '../../modules/local/unique_headers'
include { SEQKIT_TRANSLATE  } from '../../modules/local/seqkit_translate'
include { HMMER_HMMSEARCH   } from '../../modules/local/hmmer_hmmsearch'
include { MERGE_HMMER       } from '../../modules/local/merge_hmmer'


workflow HMMSEARCH {
    take:
    compactors_fasta
    spindles_fasta
    compactor_summary

    main:
    compactors_fasta.view()
    spindles_fasta.first().view()
    /*
    // Process: Add unique headers
    */
    UNIQUE_HEADERS(
        compactors_fasta.first(),
        spindles_fasta.first()
    )

    // /*
    // // Process: Translate to amino acids
    // */
    // SEQKIT_TRANSLATE(
    //     UNIQUE_HEADERS.out.fasta.flatten()
    // )

    // /*
    // // Process: Run hmmsearch
    // */
    // HMMER_HMMSEARCH(
    //     SEQKIT_TRANSLATE.out.aa,
    //     params.hmmfile
    // )

    // /*
    // // Process: Merge outputs
    // */
    // MERGE_HMMER(
    //     compactors_fasta,
    //     spindles_fasta,
    //     compactor_summary,
    //     HMMER_HMMSEARCH.out.tblout.collect()
    // )

    // emit:
    // tsv = MERGE.out.tsv

}

