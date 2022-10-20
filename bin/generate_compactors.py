#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import random
import glob
from Bio import SeqIO
import argparse
import warnings
from pandas.core.common import SettingWithCopyWarning


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fastq_samplesheet",
        type=str
    )
    parser.add_argument(
        "--anchor_length",
        type=int
    )
    parser.add_argument(
        "--epsilon",
        type=float
    )
    parser.add_argument(
        "--branch_threshold",
        type=int
    )
    parser.add_argument(
        "--recursive_depth",
        type=float
    )

    args = parser.parse_args()
    return args


def mallorn_2(dataframe, anchor_length, epsilon, N, recursive_depth, step, count_dictionary, gene_name='', effect_size=0, take_majority=True, show_leaf_only=True, valid_abundance=-1, taken_valid_abundance=False, mu_lev=0, num_annotations=0, top_ann_hit='', accept_random_nucleotide=False):

    # Protect against empty input.
    if dataframe.shape[0] == 0:
        return

    # Dictionary for utility in counts and decision history calculations.
    utility_map = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    reverse_utility = {0:'A', 1:'C', 2:'G', 3:'T', 4:'N'}

    # If recursive depth exceeded, build and report a consensus, S_value tuple.
    if recursive_depth == 0:

        # Take first index of this leaf dataframe.
        build_from = dataframe.index.tolist()[0]

        # Take anchor.
        anchor_leaf = dataframe['consensus'][build_from][0:anchor_length]

        # Take compactor.
        consensus_leaf = dataframe['consensus'][build_from]
        consensus_full_leaf = dataframe['consensus'][build_from]

        # Take both the thresholded consensus (sans majority vote) and the full recursion-depth consensus.
        sans_threshold = dataframe['path'][build_from].find('X')
        consensus_majority = False

        if sans_threshold == -1:
            consensus_leaf = dataframe['consensus'][build_from]
        if sans_threshold > 0:
            consensus_majority = dataframe['consensus'][build_from]
            valid_index = sans_threshold / 2
            consensus_leaf = consensus_majority[:int(len(anchor_leaf) + valid_index)]

        # Get the total anchors at root.
        anchor_sum = count_dictionary[anchor_leaf]

        # Get total anchors at this locus.
        bin_sum = dataframe.shape[0]

        if valid_abundance < 0:
            valid_abundance = bin_sum

        # if isinstance(top_ann_hit, float):
            # top_ann_hit = 'NaN'

        # If we haven't encountered an X, then don't forget to add the by-sample vote.
        if not taken_valid_abundance:
            # Compute the by-sample vote.
            sample_counts = dict(dataframe['sample'].value_counts())
            sample_columns = pd.read_csv('sample_specificity.tsv', engine='python', sep='\t', nrows=0).columns.tolist()[2:]
            curr_consensus = dataframe['consensus'][dataframe.index.tolist()[0]]
            curr_anchor = curr_consensus[:anchor_length]
            stork = curr_anchor + '\t' + curr_consensus + '\t'
            for column in sample_columns:
                if column in sample_counts:
                    stork += str(sample_counts[column]) + '\t'
                if column not in sample_counts:
                    stork += str(0) + '\t'
            stork = stork[:-1]
            stork += '\n'
            contribution = open('sample_specificity.tsv', 'a')
            contribution.write(stork)
            contribution.close()

        # Output this leaf dataframe as a tab-separated value table.
        naming = f"{str(bin_sum/anchor_sum)[:4]}_{anchor_sum}_{bin_sum}"
        naming_valid = f"{str(valid_abundance/anchor_sum)[:4]}_anchor_sum_{valid_abundance}"
        dataframe.to_csv(anchor_leaf + '.compactor', sep='\t')

        # Write the compactor of this leaf to an aggregate file for all anchors.
        with open('compactor.fasta', 'a') as consensus_log:
            consensus_log.write('>' + anchor_leaf + '_' + naming_valid +'\n')
            consensus_log.write(consensus_leaf+'\n')

        # Write the greedy, full-leaf on this branch.
        with open('spindles.fasta', 'a') as leaf_log:
            leaf_log.write('>' + anchor_leaf + '_' + naming +'\n')
            leaf_log.write(consensus_full_leaf +'\n')

        compactor_summary = f"{anchor_leaf}\t{consensus_leaf}\t{str(bin_sum/anchor_sum)[:4]}\t{anchor_sum}\t{str(valid_abundance/anchor_sum)[0:4]}\t{valid_abundance}\t{bin_sum}\t{dataframe['path'][build_from]}\t{consensus_full_leaf}\n"
        with open('compactor_summary.tsv', 'a') as consensus_log:
            consensus_log.write(compactor_summary)
        return

    # Proceed if we have not ended our recursion:

    # Representing A C G T N counts at current index.
    counts = np.array([0,0,0,0,0])

    # Count nucleotides.
    for i in dataframe.index.tolist():
        try:
            counts[utility_map[dataframe['sequence'][i][anchor_length + step]]] += 1
        except IndexError:
            pass

    # Compute proportions of sum.
    proportions = counts / np.sum(counts)

    # Identify for which nucleotides both conditions are true.
    bools = [counts[i] >= N and proportions[i] >= epsilon for i in np.arange(5)]
    bools_except = [counts[i] >= 5 and proportions[i] >= 0.8 for i in np.arange(5)]

    for i in np.arange(5):
        bools[i] = max(bools[i], bools_except[i])

    # Case in which no nucleotides pass thresholds; propogate most abundant anchor.
    # If this happens, we put 'X' in the decision path, rather than nucleotide abundance rank.
    if np.sum(bools) == 0 or taken_valid_abundance:
        dataframe.loc[:,'consensus'] = [dataframe['consensus'][i] + reverse_utility[np.argmax(counts)] for i in dataframe.index.tolist()]
        dataframe.loc[:,'path'] = [dataframe['path'][i] + 'X-' for i in dataframe.index.tolist()]
        # If we are encountering a failed threshold check for the first time.
        if not taken_valid_abundance:
            # Compute the by-sample vote.
            sample_counts = dict(dataframe['sample'].value_counts())
            sample_columns = pd.read_csv('sample_specificity.tsv', engine='python', sep='\t', nrows=0).columns.tolist()[2:]
            curr_consensus = dataframe['consensus'][dataframe.index.tolist()[0]]
            curr_anchor = curr_consensus[:anchor_length]
            stork = curr_anchor + '\t' + curr_consensus + '\t'
            for column in sample_columns:
                if column in sample_counts:
                    stork += str(sample_counts[column]) + '\t'
                if column not in sample_counts:
                    stork += str(0) + '\t'
            stork = stork[:-1]
            stork += '\n'
            contribution = open('sample_specificity.tsv', 'a')
            contribution.write(stork)
            contribution.close()
            return mallorn_2(dataframe, anchor_length, epsilon, N, recursive_depth - 1, step + 1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, len(dataframe.index.tolist()), True, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
        return mallorn_2(dataframe, anchor_length, epsilon, N, recursive_depth - 1, step + 1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)

    # Case in which one nucleotide passes thresholds.
    if np.sum(bools) == 1 and not taken_valid_abundance:
        dataframe.loc[:,'consensus'] = [dataframe['consensus'][i] + reverse_utility[np.argmax(bools)] for i in dataframe.index.tolist()]
        dataframe.loc[:,'path'] = [dataframe['path'][i] + '1-' for i in dataframe.index.tolist()]
        return mallorn_2(dataframe, anchor_length, epsilon, N, recursive_depth - 1, step + 1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)

    # Case in which two or more nucleotides pass thresholds.
    if np.sum(bools) >= 2 and not taken_valid_abundance:

        # Establish counts arrays for the passing and non-passing nucleotides, to be used in sampling.
        sample_space, weights, all_weights = [], [], []

        # Populate the count arrays in accordance with the 2-thresholds Truth-value vector.
        for i in np.arange(5):
            if bools[i]:
                sample_space = np.append(sample_space, i)
                weights = np.append(weights, counts[i])
                all_weights = np.append(all_weights, counts[i])
            if not bools[i]:
                all_weights = np.append(all_weights, 0)

        # Get the rank of all nucleotides as represented by their mapping integers in 'utility_map' dictionary.
        arg_sorted = np.flip(np.argsort(all_weights))
        for k in dataframe.index.tolist():

            # If we exceed the read length at this index, assign a branch as a weighted sample of valid branches.
            short_sequence = anchor_length + step > len(dataframe['sequence'][k]) - 1

            if short_sequence:
                if accept_random_nucleotide:
                    try:
                        random_nt = reverse_utility[random.choices(sample_space, weights, k=1)[0]]
                        dataframe.loc[k,'consensus'] = dataframe['consensus'][k] + random_nt
                        dataframe.loc[k,'path'] = dataframe['path'][k] + str(list(arg_sorted).index(utility_map[random_nt]) + 1) + 'M'
                    except NameError:
                        pass
                if not accept_random_nucleotide:
                    dataframe = dataframe.drop([k])

            # If we have a nucleotide at this index in the read:
            if not short_sequence:

                # And if this nucleotide passes both proportion and quantity thresholds.
                if bools[utility_map[dataframe['sequence'][k][anchor_length + step]]]:

                    # Branch in accordance with the valid nucleotide.
                    random_nt = reverse_utility[random.choices(sample_space, weights, k=1)[0]]
                    dataframe.loc[k,'consensus'] = dataframe['consensus'][k] + dataframe['sequence'][k][anchor_length+step]
                    dataframe.loc[k,'path'] = dataframe['path'][k] + str(list(arg_sorted).index(utility_map[dataframe['sequence'][k][anchor_length+step]]) + 1) + '-'

                # And if this nucleotide does not pass both proportion and quantity thresholds.
                if not bools[utility_map[dataframe['sequence'][k][anchor_length + step]]]:

                    if accept_random_nucleotide:
                        # Branch in accordance with a sampled nucleotide.
                        random_nt = reverse_utility[random.choices(sample_space, weights, k=1)[0]]
                        dataframe.loc[k,'consensus'] = dataframe['consensus'][k] + random_nt
                        dataframe.loc[k,'path'] = dataframe['path'][k] + str(list(arg_sorted).index(utility_map[random_nt]) + 1) + 'R'

                    if not accept_random_nucleotide:
                        dataframe = dataframe.drop([k])

        if take_majority:

            if np.sum(bools) == 2:

                # Recurse on tables whose reads have consistent consensus values.
                uniques = dataframe['consensus'].unique()

                recursed1 = dataframe[dataframe['consensus'] == uniques[0]]
                recursed2 = dataframe[dataframe['consensus'] == uniques[1]]

                mallorn_2(recursed1, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed2, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                return

            if np.sum(bools) == 3:

                # Recurse on tables whose reads have consistent path values.
                uniques = dataframe['consensus'].unique()

                recursed1 = dataframe[dataframe['consensus'] == uniques[0]]
                recursed2 = dataframe[dataframe['consensus'] == uniques[1]]
                recursed3 = dataframe[dataframe['consensus'] == uniques[2]]

                mallorn_2(recursed1, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed2, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed3, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                return

            if np.sum(bools) == 4:

                # Recurse on tables whose reads have consistent path values.
                uniques = dataframe['consensus'].unique()

                recursed1 = dataframe[dataframe['consensus'] == uniques[0]]
                recursed2 = dataframe[dataframe['consensus'] == uniques[1]]
                recursed3 = dataframe[dataframe['consensus'] == uniques[2]]
                recursed4 = dataframe[dataframe['consensus'] == uniques[3]]

                mallorn_2(recursed1, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed2, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed3, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed4, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                return

            if np.sum(bools) == 5:

                # Recurse on tables whose reads have consistent path values.
                uniques = dataframe['consensus'].unique()

                recursed1 = dataframe[dataframe['consensus'] == uniques[0]]
                recursed2 = dataframe[dataframe['consensus'] == uniques[1]]
                recursed3 = dataframe[dataframe['consensus'] == uniques[2]]
                recursed4 = dataframe[dataframe['consensus'] == uniques[3]]
                recursed5 = dataframe[dataframe['consensus'] == uniques[5]]

                mallorn_2(recursed1, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed2, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed3, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed4, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)
                mallorn_2(recursed5, anchor_length, epsilon, N, recursive_depth-1, step+1, count_dictionary, gene_name, effect_size, take_majority, show_leaf_only, valid_abundance, taken_valid_abundance, mu_lev, num_annotations, top_ann_hit, accept_random_nucleotide)

                return


def main():

    # get args
    args = get_args()

    print("Starting aspen")
    warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

    # Initialize sample contribution file.
    fastq_samplesheet = pd.read_csv(args.fastq_samplesheet, names=['fastq'], usecols=[0])
    fastq_ids = (
        fastq_samplesheet['fastq']
        .apply(os.path.basename)
        .str.split('.').str[0]
        .tolist()
    )
    contribution = open('sample_specificity.tsv', 'w')
    stork = ''
    for i in range(len(fastq_ids)):
        stork += fastq_ids[i] + '\t'
    stork = stork[:-1]
    stork += '\n'
    contribution.write(f'anchor\tcompactor\t{stork}')
    contribution.close()

    # Initialize compactor summary file
    header=f"anchor\tcompactor_valid\tmajority_local_proportion\tanchor_abundance\tvalid_local_proportion\tvalid_abundance\tmajority_abudance\tpath\tcompactor_majority\n"
    with open('compactor_summary.tsv', 'a') as consensus_log:
        consensus_log.write(header)

    # Process each intermediary file
    intermediarys = glob.glob("*intermediary")

    for intermediary in intermediarys:
        count_dictionary = {}
        dataframe = pd.read_csv(intermediary, names=['sequence', 'sample'], header=None, sep='\t')
        dataframe['consensus'] = [intermediary[:args.anchor_length] for i in np.arange(len(dataframe['sequence']))]
        dataframe['path'] = ['' for i in np.arange(len(dataframe['sequence']))]
        count_dictionary[intermediary[:args.anchor_length]] = dataframe.shape[0]

        mallorn_2(dataframe, args.anchor_length, args.epsilon, args.branch_threshold, args.recursive_depth, 0, count_dictionary, '', 0, True, True)

    print("Finished aspen")


main()
