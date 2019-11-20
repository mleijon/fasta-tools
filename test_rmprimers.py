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
fa.wr_fasta_file(ARGS.o)
outfi = ARGS.o.split('.')[0] + '.afa'
muscle = sub.Popen(ARGS.m + ' -in ' + ARGS.o + ' -out ' + outfi + ' -quiet')
muscle.wait()
print('Done!')