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
agent = dict()
accessions = dict()
species_set = set()
# Assumes reading 'top-files' from blast_find.py
out_file = open(ARGS.b.replace('top', 'species'), 'w')
in_file.readline()
in_file.readline()
sum_of_reads = 0
for line in in_file:
    acc = line.split('|')[0].strip()
    nr_of_reads = line.split('|')[3].strip()
    agent = {acc: nr_of_reads}
    accessions[acc]['agents'] = []
    out_file.write(acc)
# Don't consider the version of accession.
    if ma.find(acc.split('.')[0].encode('UTF-8')) == -1:  # Acc. not found.
        if 'Not found!' not in species_set:
            species_set += 'Not found!'
            accessions[acc]['species'] = 'Not found!'
            accessions[acc]['sum_of_reads'] = nr_of_reads
            accessions[acc]['agents'].append(agent)
            acc_max = acc
        else:
            sum_of_reads += nr_of_reads
            accessions[acc_max]['species'] = 'Not found!'
            accessions[acc_max]['sum_of_reads'] = sum_of_reads
            accessions[acc_max]['agents'].append(agent)
    else:
        ma.seek(ma.find(acc.split('.')[0].encode('UTF-8')))
        taxid = ma.readline().decode('UTF-8').strip().split('\t')[2]
        mn.seek(mn.find(taxid.encode('UTF-8')))
        species = mn.readline().decode('UTF-8').strip().split('|')[1]
    ma.seek(0)
    mn.seek(0)
fn.close()
fa.close()
in_file.close()
out_file.close()
