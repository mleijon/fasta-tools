#!/usr/bin/python
"""Something"""


import argparse
from multiprocessing import Process, Manager
from fasta import FastaList


def thread_work(res_lst, fasta_div):
    for fasta_t in fasta_div:
        fasta_t_seq = fasta_t.split('\n')[1]
        for fasta_s in ALL_S:
            fasta_s_id = fasta_s.split('\n')[0]
            fasta_s_seq = fasta_s.split('\n')[1]
            if fasta_s_seq.find(fasta_t_seq) != -1:
                out_start = int(fasta_s_seq.find(fasta_t_seq) - (ARGS.l - len(
                    fasta_t_seq))/2)
                out_end = int(out_start + ARGS.l+1)
                if fasta_s_seq[fasta_s_seq.find(fasta_t_seq)-12:fasta_s_seq.
                               find(fasta_t_seq)] in aivcs:
                    patho_lbl = '_HP'
                else:
                    patho_lbl = '_LP'
                out_seq = fasta_s_id + patho_lbl + '\n' + fasta_s_seq[
                    out_start:out_end] + '\n'
                res_lst.append(out_seq)


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
    ALL_S = FA_S.seq_list + FA_S.rev_comp()
    FA_T = FastaList(ARGS.t)
    FA_CS = FastaList('aivcs.fa')
    fa_t_div = FA_T.divide(ARGS.p)
    aivcs = []
    for seq in FA_CS.seq_list:
        aivcs.append(seq.split('\n')[1].strip())
    FA_OUT = open('sources.fa', 'w')
    manager = Manager()
    res_lst = manager.list()
    threads = []
    for i in range(ARGS.p):
        p = Process(target=thread_work, args=(res_lst, fa_t_div[i]))
        threads.append(p)
        p.start()
    for thread in threads:
        thread.join()
    for sequence in res_lst:
        FA_OUT.write(sequence)
    FA_OUT.close()
