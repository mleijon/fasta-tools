#!/usr/bin/python3
"""Demultiplex NGS data with multiplexed PCR amplicons listed in a fasta-file
"""
from fasta import FastaList

if __name__ == "__main__":
    import argparse

    PARSER = argparse.ArgumentParser(description='Demultiplex amplicon sequen-'
                                                 'cing NGS-data (-s) based on '
                                                 'tags listed in a separate'
                                                 'fasta-file (-t'
                                                 'Assumes the usual file '
                                                 'extensions')
    PARSER.add_argument('-s', type=str, help='NGS-data fasta/q filename',
                        required=True)
    PARSER.add_argument('-t', type=str, help='taga sequence fasta filename',
                        required=True)
    ARGS = PARSER.parse_args()
    seq_fa = FastaList(ARGS.s)
    tag_fa = FastaList(ARGS.t)
    nr_seq_demult = 0
    for tag in range(tag_fa.nr_seq):
        with open(tag_fa.id_list[tag] + '.fa', 'w') as fi:
            pass
    for seq in range(seq_fa.nr_seq):
        for tag in range(tag_fa.nr_seq):
            if tag_fa.seq_list_revc()[tag].split('\n')[1] in \
                    seq_fa.seq_list[seq].split('\n')[1]:
                with open(tag_fa.id_list[tag] + '.fa', 'a') as fi:
                    fi.write(seq_fa.id_list[seq] + seq_fa.seq_list[seq])
                    nr_seq_demult += 1
    print("{} % of {} sequences demultiplexed".format(
        round(100*nr_seq_demult/seq_fa.nr_seq), seq_fa.nr_seq))
