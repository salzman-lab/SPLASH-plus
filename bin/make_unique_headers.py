#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compactors_fasta",
        type=str
    )
    parser.add_argument(
        "--spindles_fasta",
        type=str
    )
    args = parser.parse_args()
    return args


def main():
    # Parse args
    args = get_args()
    compactor_fasta = pd.read_csv(args.compactor_fasta, sep='\t', header=None)
    spindle_fasta = pd.read_csv(args.spindles_fasta, sep='\t', header=None)

    spindles_list = list(spindle_fasta.iloc[:,0])
    spindles_file = open('spindles_id.fasta', 'a')
    for i in range(len(spindles_list)):
        if spindles_list[i][0] == '>':
            spindles_file.write(spindles_list[i] + '_ID_' + str(i) + '\n')
        else:
            spindles_file.write(spindles_list[i] + '\n')
    spindles_file.close()

    compactor_list = list(compactor_fasta.iloc[:,0])
    compactor_file = open('compactor_id.fasta', 'a')
    for i in range(len(compactor_list)):
        if compactor_list[i][0] == '>':
            compactor_file.write(compactor_list[i] + '_ID_' + str(i) + '\n')
        else:
            compactor_file.write(compactor_list[i] + '\n')
    compactor_file.close()


main()
