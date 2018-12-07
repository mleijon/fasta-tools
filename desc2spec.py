#!/usr/bin/python
""" Replace descriptions with species. """

import argparse
import mmap

PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-b', type=str, help='Fasta summary table', required=True)
PARSER.add_argument('-n', type=str, help='taxid name file', required=True)
PARSER.add_argument('-a', type=str, help='acc. to taxid file', required=True)
ARGS = PARSER.parse_args()

agents = {}
list_of_agents = []

fn = open(ARGS.n, 'r+b')
mn = mmap.mmap(fn.fileno(), 0)
fa = open(ARGS.a, 'r+b')
ma = mmap.mmap(fa.fileno(), 0)
in_file = open(ARGS.b)
out_file = open(ARGS.b.replace('top', 'species'), 'w')
out_file.write('accession\tspecies\tnr-of-reads\n')
in_file.readline()
in_file.readline()
for line in in_file:
    acc = line.split('|')[0].strip()
    nr_of_reads = line.split('|')[3].strip()
    out_file.write(acc + '\t')
    acc = acc.split('.')[0].encode('UTF-8')
    ma.seek(ma.find(acc))
    taxid = ma.readline().decode('UTF-8').strip().split('\t')[2]
    mn.seek(mn.find(taxid.encode('UTF-8')))
    species = mn.readline().decode('UTF-8').strip().split('|')[1]
    if species in agents:
        species['nr_of_reads'] += nr_of_reads
    else:
        agents['accesssion'] = acc

    out_file.write(species + '\t' + nr_of_reads + '\n')
    print(species + '\t' + nr_of_reads + '\n')
    ma.seek(0)
    mn.seek(0)
fn.close()
fa.close()
in_file.close()
out_file.close()
