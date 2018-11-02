#!/usr/bin/python3
"""Something"""


import argparse
from fasta import FastaList


PARSER = argparse.ArgumentParser(description='Finds target sequences from a\
           target (-t) fasta file in a source fasta file (-s)')
PARSER.add_argument('-t', type=str, help='target fastafile', required=True)
PARSER.add_argument('-s', type=str, help='source fastafile', required=True)
PARSER.add_argument('-l', type=str, help='the length of the extracted seq',
                    default=150)
ARGS = PARSER.parse_args()
FA_S = FastaList(ARGS.s)
FA_T = FastaList(ARGS.t)
FA_OUT = open('sources.fa', 'w')
for fasta_t in FA_T.seq_list:
    fasta_t_seq = fasta_t.split('\n')[1]
    for fasta_s in FA_S.seq_list:
        fasta_s_id = fasta_s.split('\n')[0]
        fasta_s_seq = fasta_s.split('\n')[1]
        if fasta_s_seq.find(fasta_t_seq) != -1:
            out_start = int(fasta_s_seq.find(fasta_t_seq) - (ARGS.l - len(
                fasta_t_seq))/2)
            out_end = int(out_start + ARGS.l)
            FA_OUT.write(fasta_s_id + '\n')
            FA_OUT.write(fasta_s_seq[out_start:out_end] + '\n')
FA_OUT.close()
