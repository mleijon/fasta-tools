#!/usr/bin/python
"""Metastava quickfix2 - will create a two column file with sequence name and species assignment"""

import linecache as lc
import argparse
import mmap
import collections

PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-b', type=str, help='Blast hit summary table', default='/data2/metastava/results.blast')
PARSER.add_argument('-n', type=str, help='taxid name file', default='/data2/namessci.dmp')
PARSER.add_argument('-a', type=str, help='acc. to taxid file', default='/data2/acc2taxid.dmp')
PARSER.add_argument('-o', type=str, help='Output file', default='./nonsalresults.txt')
ARGS = PARSER.parse_args()
namesf = open(ARGS.n, 'r+b')
mnamesf = mmap.mmap(namesf.fileno(), 0)
taxidf = open(ARGS.a, 'r+b')
mtaxidf = mmap.mmap(taxidf.fileno(), 0)
bl_inf = open(ARGS.b)
outf = open(ARGS.o, 'w')
result = dict()
seqname = ''
species = ''

for count, line in enumerate(bl_inf):
    if line[0:6] == 'Query=':
        seqname = line.split(' ')[1].strip()
        if 'No hits found' in lc.getline(ARGS.b, count + 7):
            species = 'NOT FOUND'
        else:
            acc = lc.getline(ARGS.b, count+9).split(' ')[0].strip()[4:]
            acc_pointer = mtaxidf.find(acc.encode('UTF-8'))
            if acc_pointer == -1:
                species = 'NOT FOUND'
            else:
                print(acc)
                mtaxidf.seek(acc_pointer)
                taxid = mtaxidf.readline().decode('UTF-8').split('\t')[1].strip()
                mnamesf.seek(mnamesf.find(taxid.encode('UTF-8')))
                species = mnamesf.readline().decode('UTF-8').split('|')[1].strip()
                result[seqname] = species

    mnamesf.seek(0)
    mtaxidf.seek(0)
ordered_result = collections.OrderedDict(sorted(result.items()))
for seq, spec in ordered_result.items():
    outf.write('%s\t%s\n' % (seq, spec))
namesf.close()
taxidf.close()
bl_inf.close()
outf.close()



