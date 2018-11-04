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
ALL_S = FA_S.seq_list + FA_S.rev_comp(FA_S.seq_list)
FA_T = FastaList(ARGS.t)
FA_OUT = open('sources.fa', 'w')
count = len(FA_T.seq_list) // 100
percent = 0
incr1 = 0
incr2 = 0
seq_found = 0
for fasta_t in FA_T.seq_list:
    incr1 += 1
    if incr1 > incr2:
        incr2 += count
        print('Calculating, %d %% completed' % min(percent,100), end='\r', flush=True)
        percent += 1
    fasta_t_seq = fasta_t.split('\n')[1]
    for fasta_s in ALL_S:
        fasta_s_id = fasta_s.split('\n')[0]
        fasta_s_seq = fasta_s.split('\n')[1]
        if fasta_s_seq.find(fasta_t_seq) != -1:
            out_start = int(fasta_s_seq.find(fasta_t_seq) - (ARGS.l - len(
                fasta_t_seq))/2)
            out_end = int(out_start + ARGS.l+1)
            FA_OUT.write(fasta_s_id + '\n')
            FA_OUT.write(fasta_s_seq[out_start:out_end] + '\n')
            seq_found += 1
FA_OUT.close()
print('')
print('sequences found: ', seq_found)
print('sequences not found: ', len(FA_S.seq_list) - seq_found)
