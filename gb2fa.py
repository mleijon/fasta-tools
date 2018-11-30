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
    for keyword in ['LOCUS', 'DEFINITION', 'ACCESSION', 'VERSION', 'KEYWORDS',
                    'SOURCE', 'ORIGIN']:
        out_fi.seek(0)
        if keyword not in out_fi.read()[:10000]:
            raise Exception('Not a Genbank File')
    out_fi.seek(0)
    return out_fi


PARSER = argparse.ArgumentParser(description='converts genbank file to fasta\
                                 file')
PARSER.add_argument('-gb', type=str, help='genbank file', required=True)
PARSER.add_argument('-d', action='store_true', help='If set indicates that gb\
                                                     is a directory that are\
                                                     scanned for gz-files')
PARSER.add_argument('-m', type=str, help='merged output file')
ARGS = PARSER.parse_args()

if ARGS.d:
    if not ARGS.gb.endswith('/'):
        ARGS.gb += '/'
    with os.scandir(ARGS.gb) as dir:
        for entry in dir:
            if entry.is_file():
                try:
                    GB_IN = rdfi(ARGS.gb + entry.name)
                except Exception:
                    print('Skipping ' + entry.name + '...')
                    continue
                print('Processing ' + entry.name + '...')
                file_list = [record + '//\n' for record in GB_IN.read().
                             split('//\n')][:-1]
                if ARGS.m is None:
                    output_file = ARGS.gb + entry.name.split('.')[0] + '.fa'
                    mode = 'w'
                else:
                    output_file = ARGS.gb + ARGS.m
                    if os.path.isfile(ARGS.gb + ARGS.m):
                        mode = 'a'
                    else:
                        mode = 'w'
                with open(output_file, mode) as fafi:
                    for record in file_list:
                        gb = GbParse(record)
                        fafi.write(gb.make_fasta() + '\n')
                if entry.name.upper().endswith('.GZ'):
                    try:
                        os.remove(ARGS.gb + entry.name[:-3])
                    except IOError:
                        pass
else:
    try:
        GB_IN = rdfi(ARGS.gb)
    except Exception:
        sys.exit('input file error')
    file_list = [record + '//\n' for record in GB_IN.read().split('//\n')][:-1]
    with open(ARGS.gb.split('.')[0] + '.fa', 'w') as fafi:
        for record in file_list:
            gb = GbParse(record)
            fafi.write(gb.make_fasta() + '\n')
