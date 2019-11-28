#!/usr/bin/python3

from fasta import FastaList
import argparse
import os
import sys
import shutil
import re
import operator
import subprocess as sub


def sep(opsys):
    if opsys == 'nt':
        return '\\'
    else:
        return '/'


def create_consensus(seqlist, gene):
    count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    consensus_id = '>' + gene + '_consensus'
    consensus_sq = ''
    for nt in range(len(seqlist[0].split('\n')[1])):
        for seq in seqlist:
            if seq.split('\n')[1][nt].upper() == 'A':
                count['A'] += 1
            elif seq.split('\n')[1][nt].upper() == 'C':
                count['C'] += 1
            elif seq.split('\n')[1][nt].upper() == 'G':
                count['G'] += 1
            elif seq.split('\n')[1][nt].upper() == 'T':
                count['T'] += 1
            else:
                pass
        consensus_sq += max(count.items(), key=operator.itemgetter(1))[0]
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    consensus = consensus_id + '\n' + consensus_sq + '\n'
    return consensus


PARSER = argparse.ArgumentParser(description='Test of primerdelete')
PARSER.add_argument('-id', type=str, help='input file directory', required=True)
PARSER.add_argument('-p', type=str, help='input primer fasta'
                                         ' file', required=True)
PARSER.add_argument('-r', type=str, help='input reference sequence fasta'
                                         ' file', required=True)
PARSER.add_argument('-od', type=str, help='output file directory', required=True)
PARSER.add_argument('-m', type=str, help='muscle path', required=True)
ARGS = PARSER.parse_args()
# Some controls of input data
if not ((ARGS.id.endswith(sep(os.name)) and
         ARGS.od.endswith(sep(os.name)))):
    sys.exit('Invalid directory name. Exits.')
if os.path.isfile(ARGS.od[:-1]):
    print('{} is a file'.format(ARGS.od[:-1]))
    sys.exit('Exits')
if os.path.isdir(ARGS.od):
    shutil.rmtree(ARGS.od)
os.mkdir(ARGS.od)
refs = FastaList(ARGS.r)
for seqfile in os.listdir(ARGS.id):
    if 'log' in seqfile:
        continue
    input_name = ARGS.id + seqfile
    output_name = ARGS.od + seqfile
    fa_in = FastaList(input_name)
    while fa_in.rmprimers(ARGS.p):
        fa_in.wr_fasta_file(output_name, ARGS.p)
        fa_in = FastaList(output_name)
    for item in refs.seq_list:
        if seqfile.split('.')[0] in item.split('\n')[0]:
            with open(output_name, 'r+') as fi:
                content = fi.read()
                fi.write(item + content)
                fi.close()
    muscle_out = output_name.split('.')[0] + '.afa'
    muscle = sub.Popen(ARGS.m + ' -in ' + output_name + ' -out ' + muscle_out +
                       ' -quiet')
    muscle.wait()
    reference_seq = ''
    cropped_alignment = FastaList(muscle_out).crop_ends()
    for item in cropped_alignment:
        if [m.start() for m in re.finditer('ref', item)]:
            reference_seq = item
    if not reference_seq:
        exit('No reference sequence. Exits')
    for i in range(len(cropped_alignment)):
        insertions = [m.start(0)
                      for m in re.finditer('-', reference_seq.split('\n')[1])]
        for index in range(len(insertions)):
            cr_al_id = cropped_alignment[i].split('\n')[0]
            cr_al_seq = cropped_alignment[i].split('\n')[1]
            new_seq = cr_al_seq[:insertions[index]] + cr_al_seq[insertions
                                                                [index] + 1:]
            cropped_alignment[i] = cr_al_id + '\n' + new_seq + '\n'
            insertions = [x - 1 for x in insertions]
    consensus = create_consensus(cropped_alignment, seqfile.split('.')[0])
    for i in range(len(cropped_alignment)):
        seq = cropped_alignment[i]
        seq_id = seq.split('\n')[0]
        seqseq = seq.split('\n')[1]
        deletions = [m.start(0) for m in re.finditer('-', seqseq)]
        if deletions:
            for index in deletions:
                seqlst = list(seqseq)
                seqlst[index] = consensus.split('\n')[1][index]
                seqseq = ''.join(seqlst)
            cropped_alignment[i] = seq_id + '\n' + seqseq + '\n'
    crfi = open(muscle_out, 'w')
    for item in cropped_alignment:
        crfi.write(item)
    crfi.close()
