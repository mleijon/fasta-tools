#!/usr/bin/python

import subprocess
import time
import datetime as dt

from fasta import FastaList


def derep(filename):
    starttime = time.time()
    fl = FastaList(filename)
    init_nrofseq = fl.nr_seq
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
    final_nrseq = len(unique_strands)
    removed_nrseq = init_nrofseq - final_nrseq
    with open(ARGS.o, 'w') as infi:
        for item in unique_strands:
            infi.write(item)
    endtime = time.time()
    if ARGS.l:
        global logfilename
        with open(logfilename, 'a') as fi:
            fi.write('*****************************'
                     '******************************\n')
            fi.write('2nd round dereplication complete in {}'.format(
                dt.timedelta(seconds=(endtime - starttime))).split(
                '.')[0])
            fi.write('\nAdditionally {} sequences removed'.
                     format(removed_nrseq))


def write_logfile():
    global logfilename
    logfilename = ARGS.o.split('.')[0] + '.log'
    with open(logfilename, 'w') as fi:
        fi.write(process1.stderr)


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
    process1 = subprocess.run([vsearch_path + 'vsearch', '--derep_prefix',
                               ARGS.f, '--sizeout',
                               '--minseqlength', ARGS.m, '--output', ARGS.o],
                              capture_output=True, text=True)
    if ARGS.l:
        global logfilename
        write_logfile()
    derep(ARGS.o)
