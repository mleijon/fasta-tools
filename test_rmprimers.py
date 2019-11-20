#!/usr/bin/python3

from fasta import FastaList
import argparse
import subprocess as sub

PARSER = argparse.ArgumentParser(description='Test of primerdelete')
PARSER.add_argument('-f', type=str, help='input file', required=True)
PARSER.add_argument('-p', type=str, help='input file', required=True)
PARSER.add_argument('-o', type=str, help='input file', required=True)
PARSER.add_argument('-m', type=str, help='muscle path', required=True)
ARGS = PARSER.parse_args()
fa = FastaList(ARGS.f)
# with open(ARGS.o, 'w') as fi:
#     for item in fa.rmprimers(ARGS.p):
#         fi.write(item)
# fi.close()
fa.rmprimers(ARGS.p)
outfi = ARGS.o.split('.')[0] + '.afa'
muscle = sub.Popen(ARGS.m + ' -in ' + ARGS.o + ' -out ' + outfi + ' -quiet')
muscle.wait()
print('Done!')