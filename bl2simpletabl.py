#!/usr/bin/python
"""Metastava quickfix2 - will create a two column file with sequence name and
species assignment"""

import linecache as lc
import multiprocessing
from multiprocessing import Pool, Manager
from functools import reduce
import argparse
import collections
import json
import blast


PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-b', type=str, help='Blast hit summary table',
                    default='test.blast')
PARSER.add_argument('-s', type=int, help='Divide blastfile into this number of '
                                         'smaller files', default=5)
PARSER.add_argument('-n', type=str, help='taxid name file', default='names.dmp')
PARSER.add_argument('-a', type=str, help='acc. to taxid file',
                    default='acc2taxid.dmp')
PARSER.add_argument('-o', type=str, help='Output file', default='nonsal.txt')
ARGS = PARSER.parse_args()
bl_inf = open(ARGS.b)
outf = open(ARGS.o, 'w')
result = dict()
acc2taxid = dict()  # {'accession':'taxid',...}
taxid2name = dict()  # {'taxid':'scientific names',...}
manager = Manager()
acc2names = manager.dict()  # {'accession':'scientific names':
blfi_list = []


def process_work(file_name):
    bl_result = dict()
    with open(file_name) as bl_fi:
        for count, row in enumerate(bl_fi):
            if row[0:6] == 'Query=':
                seq_name = row.split(' ')[1].strip()
                if 'No hits found' in lc.getline(file_name, count + 7):
                    bl_result[seq_name] = 'No hits found'
                else:
                    acc = lc.getline(file_name, count+9).\
                              split(' ')[0].strip()[4:]
                    try:
                        bl_result[seq_name] = acc2names[acc]
                    except KeyError:
                        bl_result[seq_name] = 'Accession not found'
    print('{}: Done'.format(multiprocessing.current_process()))
    return bl_result


def merge_two_dicts(x, y):
    z = x.copy()  # start with x's keys and values
    z.update(y)  # modifies z with y's keys and values & returns None
    return z


if __name__ == '__main__':
    try:
        with open('acc2name.js') as fi:
            acc2names = json.load(fi)
    except FileNotFoundError:
        with open(ARGS.a) as fi:
            for line in fi:
                acc2taxid[line.split()[1].strip()] = line.split()[2].strip()
        with open(ARGS.n) as fi:
            for line in fi:
                if line.split('|')[3].strip() == 'scientific name':
                    taxid2name[line.split('|')[0].strip()] = \
                        line.split('|')[1].strip()
        for key in acc2taxid:
            try:
                acc2names[key] = taxid2name[acc2taxid[key]]
            except KeyError:
                pass
        with open('acc2name.js', 'w') as fi:
            json.dump(acc2names, fi, indent=4)
        del taxid2name
        del acc2taxid
    print('Done loading json\n')
    blfi = blast.BlastFile(ARGS.b)
    # Split the blast file into ARGS.s smaller files and returns a list of the
    # file names
    bl_tmp_list = blast.BlastFile.split(blfi, ARGS.s)
    for item in bl_tmp_list:
        blfi_list.append(item.name)
    with Pool(processes=ARGS.s) as p:
        all_results = reduce(lambda x, y: merge_two_dicts(x, y),
                             p.starmap(process_work, blfi_list))
    ordered_result = collections.OrderedDict(sorted(all_results.items()))
    for seq, spec in ordered_result.items():
        outf.write('%s\t%s\n' % (seq, spec))
    outf.close()



