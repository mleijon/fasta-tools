#!/usr/bin/python

import subprocess
from fasta import FastaList

if __name__ == "__main__":
    import argparse

    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='fasta filename', required=True)
    PARSER.add_argument('-o', type=str, help='output fasta filename',
                        required=True)
    PARSER.add_argument('-m', type=str, help='minimun sequence length',
                        default='200')
    PARSER.add_argument('-l', action='store_true',
                        help='switch for log file output')
    ARGS = PARSER.parse_args()
    vsearch_path = '/home/micke/miniconda3/bin/'
    process = subprocess.run([vsearch_path + 'vsearch', '--derep_prefix',
                              ARGS.f, '--sizeout', '--minseqlength', ARGS.m,
                              '--output', ARGS.o], capture_output=True,
                             text=True)
    if ARGS.l:
        logfilename = ARGS.o.split('.')[0] + '.log'
        with open(logfilename, 'w') as fi:
            fi.write(process.stderr)
    fl = FastaList(ARGS.o)
    unique_strands = list()
    unique_strands.append(fl.seq_list[0])
    for item in fl:
        newstr = item.split('\n')[1].casefold()
        for i, entry in enumerate(unique_strands, start=1):
            entry_strand = entry.split('\n')[1].casefold()
            if newstr in entry_strand.casefold():
                break
            elif entry_strand.casefold() in newstr:
                unique_strands.remove(entry)
                unique_strands.append(item)
                break
            elif i == len(unique_strands):
                unique_strands.append(item)
    with open('covpo_derep2.fa', 'w') as fi:
        for item in unique_strands:
            fi.write(item)

