#!/usr/bin/python

from fasta import FastaList

if __name__ == "__main__":
    import argparse
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input fastq file', required=True)
    PARSER.add_argument('-t', type=int, help='id-tag length (default = max)', default=-1)
    ARGS = PARSER.parse_args()
    print('Loading data...',)
    fl = FastaList(ARGS.f)
    print('Done!       ')
    total_nr_of_reads = len(fl.id_list)
    unique_reads = set()
    for item in fl.seq_list:
        if ARGS.t == -1:
            unique_reads.add(item.split('\n')[1])
        else:
            unique_reads.add(item.split('\n')[1][:ARGS.t])
    unique_nr_of_reads = len(unique_reads)
    print('Total number of reads: {}'.format(total_nr_of_reads))
    print('Number of unique reads: {}'.format(unique_nr_of_reads))
    print('duplication level: {}'.format(round(total_nr_of_reads/unique_nr_of_reads, 2)))
