#!/usr/bin/python3
# Removes columns containing degenerate nucleotides from alignment and writes
# a new alignment as output file.

import argparse
from fasta import FastaList

PARSER = argparse.ArgumentParser(description='Removes columns with degenerate '
                                             'nucleotides from alignment '
                                             'and writes a new alignment as '
                                             'output file')
PARSER.add_argument('-i', type=str, help='fastq input filename', required=True)
PARSER.add_argument('-o', type=str, help='fastq output filename',
                    default='out.fa')
ARGS = PARSER.parse_args()
newseq_list = FastaList(ARGS.i)
newseq_list.rm_non_agct_columns()
newseq_list.wr_fasta_file(ARGS.o)
