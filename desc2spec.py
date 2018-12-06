#!/usr/bin/python
""" Replace descriptions with species. """

import argparse
import mmap

PARSER = argparse.ArgumentParser(description='TBD')
#PARSER.add_argument('-b', type=str, help='Fasta summary table', required=True)
PARSER.add_argument('-n', type=str, help='taxid name file', required=True)
PARSER.add_argument('-a', type=str, help='acc. to taxid file', required=True)
ARGS = PARSER.parse_args()

with open(ARGS.n, 'r+b') as fb:
    mb = mmap.mmap(fb.fileno(), 0)
    with open(ARGS.a, 'r+b') as fa:
        ma = mmap.mmap(fa.fileno(), 0)
        ma.seek(ma.find(b'X52702.1'))
        print((ma.readline()).decode('UTF-8').strip().split('\t'))
