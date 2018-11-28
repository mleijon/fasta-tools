#!/usr/bin/python3

import argparse
import sys
from genbank import GbParse


PARSER = argparse.ArgumentParser(description='converts genbank file to fasta\
                                 file')

PARSER.add_argument('-gb', type=str, help='genbank file', required=True)
ARGS = PARSER.parse_args()
try:
    GB_IN = open(ARGS.gb)
except IOError:
    sys.exit('input file error')
file_list = [record + '//' for record in GB_IN.read().split('//')][:-1]
with open(ARGS.gb.split('.')[0] + '.fa', 'w') as fafi:
    for record in file_list:
        gb = GbParse(record)
        fafi.write(gb.make_fasta() + '\n')
