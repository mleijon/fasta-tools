#!/usr/bin/python3
"""Something"""


import argparse
import threading
from fasta import FastaList


PARSER = argparse.ArgumentParser(description='Finds target sequences from a\
           target (-t) fasta file in a source fasta file (-s)')
PARSER.add_argument('-p', type=int, help='nr of processor threads', default=1)
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
res_lst = []


def thread_work(fasta_div, res_lst):
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
                sequence = fasta_s_id + patho_lbl + '\n' + fasta_s_seq[
                 out_start:out_end] + '\n'
                res_lst.append(sequence)
    event.set()
    return


threads = []
for i in range(ARGS.p):
    thread = threading.Thread(target=thread_work, args=(fa_t_div[i], res_lst))
    event = threading.Event()
    threads.append(thread)
    print('Thread: ', i)
    thread.start()
event.wait()
print(len(res_lst))
for sequence in res_lst:
    FA_OUT.write(sequence)
FA_OUT.close()
