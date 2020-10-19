#!/usr/bin/python3

def count_lines(fname):
    with open(fname) as f:
        for i, _ in enumerate(f, start=1):
            pass
    return i


def count_char(fname, ch):
    with open(fname) as f:
        return f.read().count(ch)


def reduce_fastq(fname_in, fname_out, fact):
    import linecache
    import random
    nr_of_reads = count_char(fname_in, '@')
    sample_size = round(nr_of_reads/fact)
    read_lst = list(range(nr_of_reads))
    selected_reads = random.sample(read_lst, sample_size)
    with open(fname_out, 'w') as f_out:
        for read in selected_reads:
            for i in range(4):
                f_out.write(linecache.getline(fname_in, read + i))


if __name__ == "__main__":
    import argparse
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-fq', type=str, help='input kraken file',
                        required=True)
    PARSER.add_argument('-ik', type=str, help='input kraken file',
                        required=False)
    PARSER.add_argument('-id', type=str, help='input diamond file',
                        required=False)
    PARSER.add_argument('-o', type=str, help='output kraken file',
                        required=False)
    ARGS = PARSER.parse_args()
    reduce_fastq(ARGS.fq, 'test.fastq', 100)
