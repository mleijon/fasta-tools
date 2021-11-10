#!/usr/bin/python3

import argparse
taxonomy = dict()
tax_lvls = ['k', 'd', 'p', 'c', 'o', 'f', 'g', 's']
PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-f', type=str, help='input fasta file', required=True)
PARSER.add_argument('-t', type=str, help='input taxonomy file', required=True)
PARSER.add_argument('-o', type=str, help='output file', required=True)
PARSER.add_argument('-k', type=str, help='output file', required=True)
ARGS = PARSER.parse_args()
with open(ARGS.t) as f_t, open(ARGS.f) as f_f:
    taxdata = f_t.readlines()
    seqdata = f_f.read().split('>')[1:]
for item in taxdata:
    item_acc = item.split('\t')[0]
    item_tax = item.split('\t')[1]
    item_tax = item_tax.replace(',', '.')
    item_tax = item_tax.replace('__', ':')
    item_tax = item_tax.replace(';', ',')
    item_tax = item_tax.replace(', ', ',')
    item_tax = item_tax.replace('\n', ';\n')
    for level in tax_lvls:
        if level + ':,' in item_tax:
            item_tax = item_tax.replace(level + ':,', '')
    item_tax = item_tax.replace(',s:;', ';')
    item_tax = 'tax=' + item_tax
    taxonomy[item_acc] = item_tax
with open(ARGS.o, 'w') as f:
    for item in seqdata:
        item = item.replace(item.split()[0], '>' + item.split()[0] + ';' +
                            taxonomy[item.split()[0]].strip())
        if ARGS.k in item:
            f.write(item)





