#!/usr/bin/python3
"""Something"""


import argparse
from fasta import FastaList

PARSER = argparse.ArgumentParser(description='Reverse complement a DNA strand')
PARSER.add_argument('-s', type=str, help='oligonucleotide', required=True)
ARGS = PARSER.parse_args()
STRAND = FastaList(ARGS.s)
print(STRAND.rev_comp(STRAND.seq_list))
