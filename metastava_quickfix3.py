#!/usr/bin/python
"""Metastava quickfix2 - will create a two column file with sequence name and species assignment"""

import linecache as lc
from multiprocessing import Pool, Manager
from functools import reduce
import argparse
import collections
import json
import blast


PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-b', type=str, help='Blast hit summary table', default='test2.blast')
PARSER.add_argument('-s', type=int, help='Divide blastfile into this number of smaller files', default=5)
PARSER.add_argument('-n', type=str, help='taxid name file', default='names.dmp')
PARSER.add_argument('-a', type=str, help='acc. to taxid file', default='acc2taxid.dmp')
PARSER.add_argument('-o', type=str, help='Output file', default='./nonsal.txt')
ARGS = PARSER.parse_args()
bl_inf = open(ARGS.b)
outf = open(ARGS.o, 'w')
result = dict()
acc2taxid = dict()
taxid2name = dict()
acc2names = dict()
blfi_list = []


def process_work(bl_inf):
    for count, line in enumerate(bl_inf):
        if line[0:6] == 'Query=':
            seq_name = line.split(' ')[1].strip()
            if 'No hits found' in lc.getline(ARGS.b, count + 7):
                result[seq_name] = 'No hits found'
            else:
                acc = lc.getline(ARGS.b, count+9).split(' ')[0].strip()[4:]
                try:
                    result[seq_name] = acc2names_db[acc]
                except KeyError:
                    result[seq_name] = 'Accession not found'
    return result

if __name__ == '__main__':
    try:
        with open ('acc2name.js') as fi:
            acc2names = json.load(fi)
    except FileNotFoundError:
        with open(ARGS.a) as fi:
            for line in fi:
                acc2taxid[line.split()[1].strip()] = line.split()[2].strip()
        with open(ARGS.n) as fi:
            for line in fi:
                if line.split('|')[3].strip() == 'scientific name':
                    taxid2name[line.split('|')[0].strip()] = line.split('|')[1].strip()
        for key in acc2taxid:
            try:
                acc2names[key] = taxid2name[acc2taxid[key]]
            except KeyError:
                pass
        with open ('acc2name.js', 'w') as fi:
            json.dump(acc2names, fi, indent=4)
        del taxid2name
        del acc2taxid
    blfi = blast.BlastFile(ARGS.b)
    blfi_list = blfi.split(ARGS.s)
    manager = Manager()
    acc2names_db = manager.dict(acc2names)
    with Pool(processes=ARGS.s) as p:
        all_results = reduce(lambda x,y: x + y, p.starmap(process_work, blfi_list))
    ordered_result = collections.OrderedDict(sorted(all_results.items()))
    for seq, spec in ordered_result.items():
        outf.write('%s\t%s\n' % (seq, spec))
    bl_inf.close()
    outf.close()



