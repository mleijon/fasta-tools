import argparse
from fasta import FastaList

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input fasta file',
                        required=True)
    PARSER.add_argument('-d', type=str, help='input diamond file (outfmt = 102)', required=True)
    PARSER.add_argument('--end', type=int, help='end cut', required=True)
    ARGS = PARSER.parse_args()