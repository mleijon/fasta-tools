#!/usr/bin/python3
"""Demultiplex NGS data (fa or fastq - may be gzipped) with primers listed in a
separate fasta-file
"""
from fasta import FastaList


def sep(opsys):
    """Handles the different directory separators in win/linux
     TODO: check if necessary """
    if opsys == 'nt':
        return '\\'
    else:
        return '/'


if __name__ == "__main__":
    import argparse
    import shutil
    import os
    import sys

    PARSER = argparse.ArgumentParser(description='The script substitute indels'
                                                 'with the corresponding nt in'
                                                 'a reference sequence. The'
                                                 'sequences should be in a'
                                                 'fasta file with the reference'
                                                 'sequence as the first'
                                                 'sequence')
    PARSER.add_argument('-id', type=str, help='Directory for input data',
                        required=True)
    PARSER.add_argument('-od', type=str, help='Directory for output data',
                        required=True)
    ARGS = PARSER.parse_args()
    # Some control of input file/directory names and parameter values
    if not ((ARGS.id.endswith(sep(os.name)) and
             ARGS.od.endswith(sep(os.name)))):
        sys.exit('Invalid directory name. Exits.')
    if os.path.isfile(ARGS.od[:-1]):
        print('{} is a file'.format(ARGS.od[:-1]))
        sys.exit('Exits')
    if os.path.isdir(ARGS.od):
        shutil.rmtree(ARGS.od)
    os.mkdir(ARGS.od)
    for seqfile in os.listdir(ARGS.id):
        infa = FastaList(ARGS.id + seqfile)
        ref = infa.seq_list[0]
        refid = infa.seq_list[0].split('\n')[0]
        refseq = infa.seq_list[0].split('\n')[1]
        for seq in infa.seq_list[1:2]:
            seqid = seq.split('\n')[0]
            seqseq = seq.split('\n')[1]
            for nt in range(0, len(seqseq)):
                if refseq[nt] != seqseq[nt]:
                    print(seqseq[nt])