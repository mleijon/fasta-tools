#!/usr/bin/python3

import argparse
import pickle

txid_parent = dict()
txid2taxon = dict()
acc2txid = dict()
merged = dict()
deleted = set()
txid2name = dict()
patent = '0'


def bacteria(txid):
    test = txid
    while True:
        if test == '2':
            return True
        elif test in ['0', '1', '2157', '2759', '10239', '131567']:
            return False
        else:
            test = txid_parent[test]


def archaea(txid):
    test = txid
    while True:
        if test == '2157':
            return True
        elif test in ['0', '1', '2', '2759', '10239', '131567']:
            return False
        else:
            test = txid_parent[test]


def make_taxonomy(acc):
    allowed_taxons = ['superkingdom', 'phylum', 'class', 'order', 'family',
                      'genus', 'species']
    taxon_abbr = {'superkingdom': 'k', 'phylum': 'p', 'class': 'c',
                  'order': 'o',
                  'family': 'f', 'genus': 'g', 'species': 's'}
    txid = acc2txid[acc]
    name = txid2name[txid]
    taxon = txid2taxon[txid]
    taxonomy = ''
    while True:
        if taxon in allowed_taxons:
            taxonomy = taxon_abbr[taxon] + ':' + name + ',' + taxonomy
        if name == 'Bacteria':
            taxonomy = 'tax=' + taxonomy[:-1] + ';'
            return taxonomy
        else:
            txid = txid_parent[txid]
            if txid in ['0', '1', '2157', '2759', '10239', '131567']:
                exit('Taxonomy error')
            name = txid2name[txid]
            taxon = txid2taxon[txid]


PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-a', type=str, help='input acc2txid file', required=False)
PARSER.add_argument('-n', type=str, help='input nodes file', required=False)
PARSER.add_argument('-s', type=str, help='input names file', required=False)
PARSER.add_argument('-m', type=str, help='input merged file', required=False)
PARSER.add_argument('-d', type=str, help='input delnodes file', required=False)
PARSER.add_argument('-r', type=str, help='input 16S file', required=True)
PARSER.add_argument('-o', type=str, help='output 16S Sfile', required=True)
ARGS = PARSER.parse_args()
with open(ARGS.n) as f:
    for lines in f:
        txid = lines.split('|')[0].strip()
        parent_txid = lines.split('|')[1].strip()
        taxon = lines.split('|')[2].strip()
        txid_parent[txid] = parent_txid
        txid2taxon[txid] = taxon
try:
    acc2txid = pickle.load(open('acc2txid.pickle', 'rb'))
except(OSError, IOError) as e:
    try:
        with open(ARGS.m) as f:
            for lines in f:
                merged[lines.split('|')[0].strip()] = \
                    lines.split('|')[1].strip()
        with open(ARGS.d) as f:
            deleted = set([x.strip() for x in f.read().split('|')])
        with open(ARGS.a) as f:
            for lines in f:
                if lines.split()[0].strip() == 'accession':
                    continue
                txid = lines.split()[2].strip()
                if txid == patent or txid in deleted:
                    continue
                if txid in merged.keys():
                    txid = merged[txid]
                if bacteria(txid):
                    acc2txid[lines.split()[0].strip()] = txid
        pickle.dump(acc2txid, open("acc2txid.pickle", "wb"))
    except TypeError as e:
        exit('nucl_gb.accession2taxid, merged.dmp and delnodes.dmp files needed'
             '')
try:
    txid2name = pickle.load(open('txid2name.pickle', 'rb'))
except(OSError, IOError) as e:
    try:
        with open(ARGS.s) as f:
            for lines in f:
                txid = lines.split('|')[0].strip()
                name = lines.split('|')[1].strip()
                name_class = lines.split('|')[3].strip()
                if name_class == 'scientific name' and bacteria(txid):
                    txid2name[txid] = name
        pickle.dump(txid2name, open("txid2name.pickle", "wb"))
    except TypeError as e:
        exit('names.dmp file needed')
with open(ARGS.r) as f1, open(ARGS.o, 'w') as f2:
    seqs = f1.read().split('>')[1:]
    for seq in seqs:
        defline_old = seq.split('\n', 1)[0]
        dnaseq = seq.split('\n', 1)[1]
        acc = defline_old.split('.')[0]
        try:
            tax = make_taxonomy(acc)
        except:
            continue
        defline = '>' + acc + ';' + tax + '\n'
        f2.write(defline + dnaseq)

