#!/usr/bin/python3
"""Demultiplex NGS data with multiplexed PCR amplicons listed in a fasta-file
"""
from fasta import FastaList


def sep(opsys):
    if opsys == 'nt':
        return '\\'
    else:
        return '/'


if __name__ == "__main__":
    import argparse, shutil, datetime, getpass, os, sys

    PARSER = argparse.ArgumentParser(description='Removes redundat sequences'
                                                 'after manual trimming of '
                                                 'the result of ampdemult.py')
    PARSER.add_argument('-id', type=str, help='Directory for input data',
                        required=True)
    PARSER.add_argument('-od', type=str, help='Directory for output data',
                        required=True)

    ARGS = PARSER.parse_args()

    if not ((ARGS.id.endswith(sep(os.name)) and
             ARGS.od.endswith(sep(os.name)))):
        sys.exit('Invalid directory name. Exits.')
    if os.path.isfile(ARGS.od[:-1]):
        print('{} is a file'.format(ARGS.od[:-1]))
        sys.exit('Exits')
    if os.path.isdir(ARGS.od):
        shutil.rmtree(ARGS.od)
    os.mkdir(ARGS.od)
    nr_of_files = len([name for name in os.listdir(ARGS.id) if
                       os.path.isfile(ARGS.id + name)])
    file_nr = 1
    for seqfile in os.listdir(ARGS.id):
        print('\rprocessing file {}/{}'.format(file_nr, nr_of_files), end=" ")
        inp_seq = FastaList(ARGS.id + seqfile)
        print(inp_seq.seq_list[0])
        exit()

