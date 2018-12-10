#!/usr/bin/python

"""Replace descriptions with species."""

import argparse
import mmap

PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-b', type=str, help='Fasta summary table', required=True)
PARSER.add_argument('-n', type=str, help='taxid name file', required=True)
PARSER.add_argument('-a', type=str, help='acc. to taxid file', required=True)
ARGS = PARSER.parse_args()

fn = open(ARGS.n, 'r+b')
mn = mmap.mmap(fn.fileno(), 0)
fa = open(ARGS.a, 'r+b')
ma = mmap.mmap(fa.fileno(), 0)
in_file = open(ARGS.b)
agents = dict()
accessions = dict()

# Assumes reading 'top-files' from blast_find.py
out_file = open(ARGS.b.replace('top', 'species'), 'w')
out_file.write('accession\tspecies\tnr-of-reads\n')
in_file.readline()
in_file.readline()
for line in in_file:
    acc = line.split('|')[0].strip()
    nr_of_reads = line.split('|')[3].strip()
    agents[acc]['nr_of_reads'] = nr_of_reads

    out_file.write(acc)
# Don't consider the version of accession.
    if ma.find(acc.split('.')[0].encode('UTF-8')) == -1:  # Acc. not found.
        out_file.write(' not found!\n')
        print(acc + ' not found!\n')
    else:
        ma.seek(ma.find(acc.split('.')[0].encode('UTF-8')))
        taxid = ma.readline().decode('UTF-8').strip().split('\t')[2]
        if mn.find(taxid.encode('UTF-8')) == -1:
            out_file.write('taxid: ' + taxid + ' not found!\n')
            print('taxid: ' + taxid + ' not found!\n')
        else:
            mn.seek(mn.find(taxid.encode('UTF-8')))
            species = mn.readline().decode('UTF-8').strip().split('|')[1]
            out_file.write(species + '\t' + nr_of_reads + '\n')
            print(species + '\t' + nr_of_reads + '\n')
    ma.seek(0)
    mn.seek(0)
fn.close()
fa.close()
in_file.close()
out_file.close()
