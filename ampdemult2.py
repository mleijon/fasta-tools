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
    import copy
    import collections
    import argparse, shutil, os, sys

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
    for seqfile in os.listdir(ARGS.id):
        sample = dict()
        sample_reduced = dict(dict())
        inp_seqs = FastaList(ARGS.id + seqfile)
        for seq in inp_seqs.seq_list:
            seq_name = seq.split('_count:')[0][1:]
            seq_count = int(seq.split('_count:')[1].split('_')[0])
            dna_seq = seq.split('\n')[1]
            sample_name = seq_name.rsplit('_', 1)[0]
            sample[seq_name] = [sample_name, dna_seq, seq_count]
        for key1 in sample:
            found = False
            if sample[key1][0] in sample_reduced.keys():
                for key2 in sample_reduced[sample[key1][0]]:
                    if sample[key1][1] ==\
                            sample_reduced[sample[key1][0]][key2][1]:
                        sample_reduced[sample[key1][0]][key2][2] +=\
                            sample[key1][2]
                        found = True
                if not found:
                    sample_reduced[sample[key1][0]][key1] = sample[key1]
            else:
                sample_reduced[sample[key1][0]] = {key1: sample[key1]}
        sample_tmp = dict()
        for key in sample_reduced:
            n = key.split('_')
            name = n[0] + '_' + n[1].zfill(2) + '_' + n[2].\
                zfill(2) + '_' + n[3].zfill(2)
            sample_tmp[name] = sample_reduced[key]
        sample_reduced = copy.deepcopy(sample_tmp)
        samples_red_ord = collections.OrderedDict(sorted(
            sample_reduced.items()))
        with open(ARGS.od + seqfile, 'w') as fi:
            for key1 in samples_red_ord:
                for key2 in samples_red_ord[key1]:
                    n = key2.split('_')
                    name = n[0] + '.' + n[1].zfill(2) + '.' + n[2].\
                        zfill(2) + '.' + n[3].zfill(2) + '_' + n[4][n[4].rfind('0') + 1:]
                    dna_seq = samples_red_ord[key1][key2][1]
                    seq_count = samples_red_ord[key1][key2][2]
                    fi.write('>{}_c:{}\n{}\n'.
                             format(name, seq_count, dna_seq))




