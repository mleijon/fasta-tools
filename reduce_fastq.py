#!/usr/bin/python3


def count_char(fname, ch):
    with open(fname) as f:
        return f.read().count(ch)


def reduce_fastq(fname_in, fact):
    import linecache
    import random
    nr_of_reads = count_char(fname_in, '@')
    print('{} sequence reads'.format(nr_of_reads))
    sample_size = round(nr_of_reads/fact)
    read_lst = list(range(1, nr_of_reads, 4))
    selected_reads = random.sample(read_lst, sample_size)
    r1_new = fname_in.rpartition('/')[0] + fname_in.rpartition('/')[1] + 'r_' + fname_in.rpartition('/')[2]
    r2_new = r1_new.replace('_R1_', '_R2_')
    with open(r1_new, 'w') as r1_out, open(r2_new, 'w') as r2_out:
        for read in selected_reads:
            print('reads: {}'.format(read))
            for i in range(4):
                r1_out.write(linecache.getline(fname_in, read + i))
                linecache.clearcache()
                r2_out.write(linecache.getline(fname_in.replace('_R1_', '_R2_'), read + i))
                linecache.clearcache()


if __name__ == "__main__":
    import argparse
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input R1 fastq file', required=True)
    ARGS = PARSER.parse_args()
    reduce_fastq(ARGS.f, 1000)
