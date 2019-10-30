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
    sample = dict()
    sample_reduced = dict()
    for seqfile in os.listdir(ARGS.id):
        print('\rprocessing file {}/{}'.format(file_nr, nr_of_files), end=" ")
        inp_seqs = FastaList(ARGS.id + seqfile)
        for seq in inp_seqs.seq_list:
            seq_name = seq.split('_count:')[0][1:]
            seq_count = int(seq.split('_count:')[1].split('_')[0])
            seq_length = int(seq.split(':')[2].split('_')[0])
            dna_seq = seq.split('\n')[1]
            sample_name = seq_name.rsplit('_', 1)[0]
            sample[seq_name] = (sample_name, dna_seq, seq_count, seq_length)
        file_nr += 1


