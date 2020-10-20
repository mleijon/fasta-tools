#!/usr/bin/python3

def count_lines(fname):
    with open(fname) as f:
        for i, _ in enumerate(f, start=1):
            pass
    return i


def count_char(fname, ch):
    with open(fname) as f:
        return f.read().count(ch)


if __name__ == "__main__":
    import argparse
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-fq', type=str, help='input fastq file',
                        required=True)
    PARSER.add_argument('-ik', type=str, help='input kraken file',
                        required=False)
    PARSER.add_argument('-id', type=str, help='input diamond file',
                        required=False)
    PARSER.add_argument('-o', type=str, help='output kraken file',
                        required=False)
    ARGS = PARSER.parse_args()
    print(count_char(ARGS.fq, '@'))

