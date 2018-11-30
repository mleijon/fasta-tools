#!/usr/local/miniconda3/bin/python

import argparse as ap
from fasta import FastaList

PARSER = ap.ArgumentParser(description='Check fastafile')
PARSER.add_argument('-f', type=str, help='fasta file', required=True)
ARGS = PARSER.parse_args()

fl = FastaList(ARGS.f)
for seq in fl.seq_list:
    for nt in seq.split('\n')[1]:
        if nt not in ['A', 'C', 'G', 'T', 'R', 'W', 'K', 'Y', 'S', 'N', 'M',
                      'V', 'D', 'H', 'B']:
            #print(nt)
            print(seq.split('\n')[0])
