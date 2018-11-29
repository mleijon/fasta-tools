#!/usr/local/miniconda3/bin/python

import argparse
import sys
import os
import gzip
from genbank import GbParse


def rdfi(in_fi):
    """ Reads and return a genbank file unzipped if necessary"""
    if in_fi.upper().endswith('GZ'):
        out_fi = open(in_fi[:-3], 'w')
        with gzip.open(in_fi, 'rt') as fi:
            out_fi.write(fi.read())
        out_fi = open(in_fi[:-3])
    else:
        out_fi = open(in_fi)
    return out_fi


PARSER = argparse.ArgumentParser(description='converts genbank file to fasta\
                                 file')
PARSER.add_argument('-gb', type=str, help='genbank file', required=True)
PARSER.add_argument('-d', action='store_true', help='If set indicates that gb\
                                                     is a directory that are\
                                                     scanned for gz-files')
ARGS = PARSER.parse_args()

if ARGS.d:
    with os.scandir(ARGS.gb) as dir:
        for entry in dir:
            if entry.name.endswith('.gz') and entry.is_file():
                try:
                    GB_IN = rdfi(ARGS.gb + entry.name)
                    print('Processing ' + entry.name + '...', end='\r',
                          flush=True)
                except IOError:
                    sys.exit('input file error')
                file_list = [record + '//\n' for record in GB_IN.read().
                             split('//\n')][:-1]
                with open(ARGS.gb + entry.name.split('.')[0] + '.fa', 'w') as\
                        fafi:
                    for record in file_list:
                        gb = GbParse(record)
                        fafi.write(gb.make_fasta() + '\n')
else:
    try:
        GB_IN = rdfi(ARGS.gb)
    except IOError:
        sys.exit('input file error')
    file_list = [record + '//\n' for record in GB_IN.read().split('//\n')][:-1]
    with open(ARGS.gb.split('.')[0] + '.fa', 'w') as fafi:
        for record in file_list:
            gb = GbParse(record)
            fafi.write(gb.make_fasta() + '\n')
