#!/usr/bin/python3
"""Something"""


import argparse
from fasta import FastaList


def find_target():
    pass


PARSER = argparse.ArgumentParser(description='Finds target sequences from a\
           target (-t) fasta file in a source fasta file (-s)')
PARSER.add_argument('-t', type=str, help='target fastafile', required=True)
PARSER.add_argument('-s', type=str, help='source fastafile', required=True)
PARSER.add_argument('-l', type=str, help='the length of the extracted seq',
                    default=150)
ARGS = PARSER.parse_args()
FA_S = FastaList(ARGS.s)
FA_S_seqs = list(map(lambda x: x.replace('\n', ''), FA_S.seq_list))
FA_T = FastaList(ARGS.t)
FA_T_seqs = list(map(lambda x: x.replace('\n', ''), FA_T.seq_list))
FA_OUT = open('sources.fa', 'w')
seq_len = ARGS.l
for seqid_t in FA_T.id_list:
    for seqid_s in FA_S.id_list:
        pass
