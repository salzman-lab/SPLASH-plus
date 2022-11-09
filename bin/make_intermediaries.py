#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import gzip
import random
import time
import sys
import argparse
import logging
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fastq",
        type=str
    )
    parser.add_argument(
        "--fastq_id",
        type=str
    )
    parser.add_argument(
        "--anchors_file",
        type=str
    )
    parser.add_argument(
        "--anchor_length",
        type=int
    )
    parser.add_argument(
        "--num_reads",
        type=int
    )
    parser.add_argument(
        "--max_anchor_reads",
        type=int
    )
    parser.add_argument(
        "--outfile_log",
        type=str
    )
    args = parser.parse_args()
    return args


def make_intermediaries(anchor_list, args):

    # Construct anchor dictionaries.
    zeroes = [0 for i in range(len(anchor_list))]
    count_dictionary = dict(zip(anchor_list, zeroes))
    empty_strings = [[] for i in range(len(anchor_list))]
    read_dictionary = dict(zip(anchor_list, empty_strings))

    # Set threshold for reads per fastq: default to 100 million / number of samples.
    num_reads = args.num_reads

    # Utility.
    num_anchors = len(anchor_list)

    # Create variables for performance tracking.
    local_reads = 0
    written_reads = 0

    # Log the start time.
    logging.info("Starting to make intermediaries.")
    start_time = time.time()

    # Create dictionary for per-sample anchor count thresholding.
    fastq_specific_count = dict(zip(anchor_list, zeroes))

    # Open and access records using BioPython.
    with gzip.open(args.fastq, mode="rt") as handle:
        for record in SeqIO.parse(handle, 'fastq'):

            # Break if we've exceeded the number of reads we wish to visit.
            # Break if we've fully populated each of our anchor-read lists.
            if local_reads >= num_reads or written_reads == num_anchors * args.max_anchor_reads:
                break

            # If we haven't exceeded our read bounds, proceed and increment counter.
            if not local_reads >= num_reads and not written_reads == num_anchors * args.max_anchor_reads:
                local_reads += 1
                read = str(record.seq)

                # Index only to the read_length - 50th location in the read.
                for i in np.arange(len(read) - 50):

                    # Establish anchor so as to avoid slicing it from read repeatedly.
                    anchor = read[i:i+args.anchor_length]

                    # Check if we've exceeded count for an anchor.
                    # This simultaneously allows us to pass a read window if it's not an anchor.
                    if anchor in fastq_specific_count:

                        # Increment the count of this anchor in this sample.
                        fastq_specific_count[anchor] = fastq_specific_count[anchor] + 1

                        # If we are at or below our count limit:
                        if fastq_specific_count[anchor] <= args.max_anchor_reads:

                            # Add read to dictionary and increment sample-specific written read count.
                            read_dictionary[anchor].append(f"{read[i:]}\t{args.fastq_id}\n")
                            written_reads += 1

        # Log filename, reads, writes, and execution time.
        logging.info("Finished making intermediaries.")
        end_time = time.time() - start_time
        logging.info(f"{args.fastq}\t{local_reads}\t{written_reads}\t{end_time}")

    logging.info("Starting to write out.")
    # Write the resulting read lists to the appropriate intermediate files.
    for anchor in read_dictionary.keys():
        with open(f"{anchor}_{args.fastq_id}.intermediary", 'a') as anchor_log:
            for item in read_dictionary[anchor]:
                anchor_log.write(item)
        logging.info(f"\tDone with {anchor}.")
    logging.info("Finished writing out.")

    return


def main():
    # Parse args
    args = get_args()

    # Set up logging
    logging.basicConfig(
        filename=args.outfile_log,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%H:%M:%S',
        level=logging.DEBUG
    )

    # Parse anchor list
    anchor_list = (
        pd.read_csv(args.anchors_file, sep='\t', usecols=[0], names=['anchor'])
        .sort_values(['anchor'])
        ['anchor']
        .tolist()
    )

    # Make intermediate files
    make_intermediaries(anchor_list, args)


main()
