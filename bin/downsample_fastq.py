#!/usr/bin/env python3

import argparse
import random
import gzip
from Bio import SeqIO

def downsample_fastq(input_file, downsample_rate, output_file):
    """
    Downsample a single FASTQ file.

    Parameters:
    input_file (str): Path to the input FASTQ file.
    downsample_rate (float): Downsample rate (between 0 and 1).
    output_file (str): Path to the output FASTQ file.

    Returns:
    None
    """
    is_gzipped = input_file.endswith(".gz")
    open_func = gzip.open if is_gzipped else open

    with open_func(input_file, "rt") as input_handle, open_func(output_file, "wt") as output_handle:
        for record in SeqIO.parse(input_handle, "fastq"):
            if random.random() < downsample_rate:
                SeqIO.write(record, output_handle, "fastq")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Downsample a single FASTQ file.")
    parser.add_argument("input_file", help="Path to the input FASTQ file.")
    parser.add_argument("downsample_rate", type=float, help="Downsample rate (between 0 and 1).")
    parser.add_argument("output_file", help="Path to the output FASTQ file.")
    args = parser.parse_args()

    downsample_fastq(args.input_file, args.downsample_rate, args.output_file)