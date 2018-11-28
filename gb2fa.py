#!/usr/bin/python3

import argparse
import sys
import gzip
from genbank import GbParse

def rdfi(in_fi):
    """ Reads and return a genbank file unzipped if necessary"""
    if in_fi.upper().endswith('GZ'):
        out_fi = open(in_fi[:-3],'w')
        with gzip.open(in_fi, 'rt') as fi:
            out_fi.write(fi.read())
        out_fi = open(in_fi[:-3])
    else:
        out_fi = open(in_fi)
    return out_fi

PARSER = argparse.ArgumentParser(description='converts genbank file to fasta\
                                 file')
PARSER.add_argument('-gb', type=str, help='genbank file', required=True)
ARGS = PARSER.parse_args()
try:
    GB_IN = rdfi(ARGS.gb)
except IOError:
    sys.exit('input file error')
file_list = [record + '//' for record in GB_IN.read().split('//')][:-1]
with open(ARGS.gb.split('.')[0] + '.fa', 'w') as fafi:
    for record in file_list:
        gb = GbParse(record)
        fafi.write(gb.make_fasta() + '\n')
