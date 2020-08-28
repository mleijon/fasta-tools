#!/usr/bin/python3

from collections import Counter
import argparse

ROOT = '1'
CELLULAR_ORGANISM = '131567'
ARCHAEA = '2157'
BACTERIA = '2'
EUCARYOTA = '2759'
VIRUSES = '10239'
UNCLASSIFIED = '12908'
OTHERSEQ = '28384'


def count_rows(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def count_0and1(fname):
    count0 = 0
    count1 = 0
    with open(fname) as f:
        for row in f:
            item = row.split()
            if item[1].strip() == '0':
                count0 += 1
            elif item[1].strip() == '1':
                count1 += 1
    return count0, count1


def count_txclass(fname, taxclass):
    if taxclass == ARCHAEA:
        taxclass_set = archaea_taxids
    elif taxclass == BACTERIA:
        taxclass_set = bacteria_taxids
    elif taxclass == EUCARYOTA:
        taxclass_set = eucaryota_taxids
    elif taxclass == VIRUSES:
        taxclass_set = viruses_taxids
    elif taxclass == UNCLASSIFIED:
        taxclass_set = unclassified_taxids
    elif taxclass == OTHERSEQ:
        taxclass_set = otherseq_taxids
    else:
        taxclass_set = 'NONE'
        print('No such taxonomy class. Exits.')
        exit()
    count = 0
    with open(fname) as f:
        for row in f:
            taxid = row.split()[1].strip()
            if taxid in taxclass_set:
                count += 1
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxclass_set):
                count += 1
    return count


def sort_taxids(fname, taxclass):
    taxclass_set = set()
    with open(fname) as f:
        for line in f:
            line = line.replace('|', ' ')
            item = list(map(str.strip, line.split()[0:3]))
            if taxclass in item:
                taxclass_set.add(item[0])
    return taxclass_set


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input diamond file',
                        required=True)
    PARSER.add_argument('-l', type=str, help='output log file',
                        required=False, default='daa2species.log')
    ARGS = PARSER.parse_args()

sciNames = dict()
mergedTaxids = dict()
deletedTaxids = set()
Results = []
nr_of_reads = count_rows(ARGS.f)
nr_of_nohits = count_0and1(ARGS.f)
archaea_taxids = sort_taxids("taxidlineage.dmp", ARCHAEA)
bacteria_taxids = sort_taxids("taxidlineage.dmp", BACTERIA)
eucaryota_taxids = sort_taxids("taxidlineage.dmp", EUCARYOTA)
viruses_taxids = sort_taxids("taxidlineage.dmp", VIRUSES)
unclassified_taxids = sort_taxids("taxidlineage.dmp", UNCLASSIFIED)
otherseq_taxids = sort_taxids("taxidlineage.dmp", OTHERSEQ)

names_file = open("names.dmp")
merge_file = open("merged.dmp")
delete_file = open("delnodes.dmp")
lineage_file = open("taxidlineage.dmp")

# Creates dictionary with merged taxids
for line in merge_file:
    item = line.split('|')
    mergedTaxids.update({item[0].strip(): item[1].strip()})

# Creates a set with deleted taxids
for line in delete_file:
    item = line.split('|')[0].strip()
    deletedTaxids.add(item)

merge_file.close()
delete_file.close()

# Creates a dictionary with key:values equal taxid:scientific name
for line in names_file:
    item = line.split('|')
    if "scientific" in item[3]:
        sciNames.update({item[0].strip(): item[1].strip()})
names_file.close()

with open(ARGS.f) as f:
    deleted_count = 0
    cellorg_count = 0
    for line in f:
        item = line.split()[1].strip()
        if item != '0':
            if item in mergedTaxids.keys():
                Results += [sciNames[mergedTaxids[item]]]
            elif item in deletedTaxids:
                deleted_count += 1
                Results += [item + ' :deleted']
            elif item == CELLULAR_ORGANISM:
                cellorg_count += 1
            else:
                Results += [sciNames[item]]

Results = Counter(Results).most_common()
SUM_IN_DB = 0
SUM_IN_DB += count_txclass(ARGS.f, ARCHAEA)
SUM_IN_DB += count_txclass(ARGS.f, BACTERIA)
SUM_IN_DB += count_txclass(ARGS.f, EUCARYOTA)
SUM_IN_DB += count_txclass(ARGS.f, VIRUSES)
SUM_IN_DB += count_txclass(ARGS.f, UNCLASSIFIED)
SUM_IN_DB += count_txclass(ARGS.f, OTHERSEQ)
SUM_IN_DB += cellorg_count

print('Total nr of reads: {:>31}'.format(nr_of_reads))
print('\nNo hits to db: {:>35} ({} %)'.format(nr_of_nohits[0],
                                            (round(100 * (nr_of_nohits[0] /
                                                          nr_of_reads), 2))))
print('Could not be classified by db: {:>19} ({} %)'.format(
    nr_of_nohits[1], (round(100 * (nr_of_nohits[1] / nr_of_reads), 2))))
print('Taxid deleted from db: {:>27} ({} %)'.format(
    deleted_count, (round(100 * (deleted_count / nr_of_reads), 2))))
print('\nClassified by db: {:>32} ({} %)'.format(
    SUM_IN_DB, (round(100 * (SUM_IN_DB / nr_of_reads), 2))))
print('------------------------------------------------------------')
print('CELLULAR ORGANISM: {:>31} ({} %)'.format(cellorg_count, (
    round(100 * (cellorg_count / nr_of_reads), 2))))
print('ARCHAEA: {:>41} ({} %)'.format(count_txclass(ARGS.f, ARCHAEA), (
    round(100 * (count_txclass(ARGS.f, ARCHAEA) / nr_of_reads), 2))))
print('BACTERIA: {:>40} ({} %)'.format(count_txclass(ARGS.f, BACTERIA), (
    round(100 * (count_txclass(ARGS.f, BACTERIA) / nr_of_reads), 2))))
print('EUCARYOTA: {:>39} ({} %)'.format(count_txclass(ARGS.f, EUCARYOTA), (
    round(100 * (count_txclass(ARGS.f, EUCARYOTA) / nr_of_reads), 2))))
print('VIRUS: {:>43} ({} %)'.format(count_txclass(ARGS.f, VIRUSES), (
    round(100 * (count_txclass(ARGS.f, VIRUSES) / nr_of_reads), 2))))
print('UNCLASSIFIED: {:>36} ({} %)'.format(count_txclass(ARGS.f,
                                                         UNCLASSIFIED), (round(
    100 * (count_txclass(ARGS.f, UNCLASSIFIED)
           / nr_of_reads), 2))))
print('OTHERSEQ: {:>40} ({} %)'.format(count_txclass(ARGS.f, OTHERSEQ), (
    round(100 * (count_txclass(ARGS.f, OTHERSEQ) / nr_of_reads), 2))))
