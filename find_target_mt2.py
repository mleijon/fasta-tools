#!/usr/bin/python3
"""Something"""


import argparse
import json
from multiprocessing import Pool, Manager
from functools import reduce
from fasta import FastaList


def process_work(fasta_div):
    """Creates AIV HA gene fragments around the CS"""
    global ALL_S
    tmp_lst = []
    counter = 0
    found = set()
    for fasta_t in fasta_div:
        fasta_t_seq = fasta_t.split('\n')[1]
        ALL_S = [seq for seq in ALL_S if seq not in found]
        found.clear()
        for fasta_s in ALL_S:
            if fasta_s[1].find(fasta_t_seq) != -1:
                out_start = int(fasta_s[1].find(fasta_t_seq) - (ARGS.l - len(
                    fasta_t_seq))/2)
                out_end = int(out_start + ARGS.l+1)
                if fasta_s[1][fasta_s[1].find(fasta_t[1])-12:fasta_s[1].
                              find(fasta_t[1])] in aivcs:
                    patho_lbl = '_HP'
                else:
                    patho_lbl = '_LP'
                out_seq = fasta_s[0] + patho_lbl + '\n' + fasta_s[1][
                    out_start:out_end] + '\n'
                tmp_lst.append(out_seq)
                found.add(fasta_s)
        counter += 1
        # print('{:03} % completed'.format(int((100*counter/len(fasta_div)))),
        #       end='\r', flush=True)
    return tmp_lst


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Finds target sequences from\
    a target (-t) fasta file in a source fasta file (-s)')
    PARSER.add_argument('-p', type=int, help='nr of processor threads',
                        default=1)
    PARSER.add_argument('-t', type=str, help='target fastafile', required=True)
    PARSER.add_argument('-s', type=str, help='source fastafile', required=True)
    PARSER.add_argument('-l', type=int, help='the length of the extracted seq',
                        default=150)
    ARGS = PARSER.parse_args()
    FA_S = FastaList(ARGS.s)
    manager = Manager()
    ALL_S = manager.list([(fasta_s.split('\n')[0], fasta_s.split('\n')[1])
                          for fasta_s in FA_S.seq_list])
    FA_T = FastaList(ARGS.t)
    fa_t_div = FA_T.divide(ARGS.p)
    nr_s_seq = len(ALL_S)
    nr_t_seq = len(FA_T.seq_list)
    print('searching {} subsequences in {} sequences'.format(
     nr_t_seq, nr_s_seq))
    argument = [(x, y) for x in fa_t_div for y in [ALL_S]]
    aivcs = []
    re_lst = []
    try:
        with open('cs.js') as cs_file:
            aivcs = json.load(cs_file)
    except FileNotFoundError:
        FA_CS = FastaList('aivcs.fa')
        for seq in FA_CS.seq_list:
            aivcs.append(seq.split('\n')[1].strip())
        with open('cs.js', 'w') as cs_file:
            json.dump(aivcs, cs_file, indent=4)
    FA_OUT = open('sources.fa', 'w')
    with Pool(processes=ARGS.p) as p:
        re_lst = reduce(lambda x, y: x + y, p.map(process_work, fa_t_div))
    p.close()
    print('\n\n\n{:5} sequence(s) found'.format(len(re_lst)))
    print('{:5} sequence(s) not found'.format(
     len(FA_S.seq_list) - len(re_lst)))
    for sequence in re_lst:
        FA_OUT.write(sequence)
    FA_OUT.close()
