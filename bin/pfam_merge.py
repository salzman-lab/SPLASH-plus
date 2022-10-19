#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compactor_summary",
        type=str
    )
    parser.add_argument(
        "--compactors_fasta",
        type=str
    )
    parser.add_argument(
        "--spindles_fasta",
        type=str
    )
    parser.add_argument(
        "--compactors_tblout",
        type=str
    )
    parser.add_argument(
        "--spindles_tblout",
        type=str
    )
    parser.add_argument(
        "--outfile",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    # Parse args
    args = get_args()

    # Load appropriate files from compactor and Pfam output.

    compactor_summary = pd.read_csv(args.compactor_summary, sep='\t')
    compactor_fasta = pd.read_csv(args.compactors_fasta, sep='\t', header=None)
    spindle_fasta = pd.read_csv(args.compactor_summary, sep='\t', header=None)
    compactor_pfam = pd.read_csv(args.compactors_tblout, engine='c', sep='\t')
    spindle_pfam = pd.read_csv(args.spindles_tblout, engine='c', sep='\t')


    # Hack Pfam output for spindles and compactors into a table.

    fields = ['target_name', 'accession', 'query_name', 'accession2', 'full_seq_evalue', 'full_seq_score', 'full_seq_bias', 'best_1_domain_evalue', 'best_1_domain_score', 'best_1_domain_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']
    fields_spindle = [i + '_spindle' for i in fields if i != 'target_name']
    fields_spindle.insert(0,'target_name')

    pfam_list_compactor = [i for i in list(compactor_pfam.iloc[:,0]) if str(i)[0] in 'ACGT']
    pfam_list_spindles =  [i for i in list(spindle_pfam.iloc[:,0]) if str(i)[0] in 'ACGT']

    compactor_ok = [ i.split(' ') for i in pfam_list_compactor]
    compactor_okok = []
    for i in compactor_ok:
        compactor_okok.append([j for j in i if j])
    compactor_pfam_structured = pd.DataFrame(compactor_okok, columns = fields)
    compactor_pfam_structured.loc[:,"target_name"] = [i.split("_frame")[0] for i in list(compactor_pfam_structured.loc[:,"target_name"])]

    spindle_ok = [ i.split(' ') for i in pfam_list_spindles]
    spindle_okok = []
    for i in spindle_ok:
        spindle_okok.append([j for j in i if j])
    spindle_pfam_structured = pd.DataFrame(spindle_okok, columns = fields_spindle)
    spindle_pfam_structured.loc[:,"target_name"] = [i.split("_frame")[0] for i in list(spindle_pfam_structured.loc[:,"target_name"])]


    # Generate a marker (FASTA header) to guide Pfam table in merge with summary file.

    compactors = list(compactor_fasta.iloc[:,0])
    comp_headers = [i[1:] for i in compactors if '>' in i]
    comp_seqs = [i for i in compactors if i[0] in 'ACGT']
    compactor_fasta_df = pd.DataFrame({'target_name' : comp_headers, 'compactor_valid' : comp_seqs})

    summary_fasta_merged = compactor_summary.merge(compactor_fasta_df, how='left', on='compactor_valid')


    # Merge structured Pfam results into the appropriate files.

    compactor_summary_pfam_merged = summary_fasta_merged.merge(compactor_pfam_structured, how='left', on='target_name')
    spindle_summary_pfam_merged = compactor_summary_pfam_merged.merge(spindle_pfam_structured, how='left', on='target_name')

    # Write to TSV.

    spindle_summary_pfam_merged.to_csv('compactor_summary_pfam_merged.tsv', sep='\t')


main()
