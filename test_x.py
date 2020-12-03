#!/usr/bin/python

import subprocess
import time
import datetime as dt

from fasta import FastaList


def derep():
    starttime = time.time()
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


def rmdegenerates():
    pass


def write_logfile():
    with open(logfilename, 'w') as fi:
        fi.write(proc1.stderr.read())
        fi.write(proc2.stderr.read())


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

    fl = FastaList(ARGS.f)

    vsearch_path = '/home/micke/miniconda3/bin/'
    proc1args = [vsearch_path + 'vsearch --derep_prefix ' + ARGS.f +
                 ' --sizeout --minseqlength ' + ARGS.m + ' --output '
                 + '/dev/stdout']
    proc2args = [vsearch_path + 'vsearch', '--sortbylength', '-', '--output',
                 ARGS.o]
    proc1 = subprocess.Popen(proc1args, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True)
    proc2 = subprocess.Popen(proc2args, stdin=proc1.stdout,
                             stderr=subprocess.PIPE, text=True)
    proc2.wait()
    if ARGS.l:
        logfilename = ARGS.o.split('.')[0] + '.log'
        write_logfile()
    derep()
