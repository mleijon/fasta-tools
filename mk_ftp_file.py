#!/usr/bin/python3

import argparse

PARSER = argparse.ArgumentParser(description='output filename')
PARSER.add_argument('-f', type=str, help='output filename', required=True)
PARSER.add_argument('-c', type=int, help='nr of files', required=True)
ARGS = PARSER.parse_args()
startstring = "ftp://ftp.ncbi.nlm.nih.gov/ncbi-asn1/gbvrlX.aso.gz\n"
with open(ARGS.f, 'w') as fi:
    for i in range(1, ARGS.c + 1):
        fi.write(startstring.replace('X', str(i)))
fi.close()
