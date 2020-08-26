#!/usr/bin/python3

from collections import Counter
import argparse


def count_rows(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input diamond file', required=True)
    PARSER.add_argument('-l', type=str, help='output log file', required=False, default='daa2species.log')
    ARGS = PARSER.parse_args()


sciNames = dict()
mergedTaxids = dict()
deletedTaxids = set()
Results = []
nr_of_reads = count_rows(ARGS.f)
print(nr_of_reads)
exit()
names_file = open("names.dmp")
merge_file = open("merged.dmp")
delete_file = open("delnodes.dmp")
lineage_file = open("taxidlineage.dmp")
for line in merge_file:
    item = line.split('|')
    mergedTaxids.update({item[0].strip(): item[1].strip()})
for line in delete_file:
    item = line.split('|')[0].strip()
    deletedTaxids.add(item)
merge_file.close()
daa_file = open(ARGS.f)
for line in names_file:
    item = line.split('|')
    if "scientific" in item[3]:
        sciNames.update({item[0].strip(): item[1].strip()})
names_file.close()
for line in daa_file:
    item = line.split()[1].strip()
    if item != '0':
        if item in mergedTaxids.keys():
            Results += [sciNames[mergedTaxids[item]]]
        elif item in deletedTaxids:
            Results += [item + ' :deleted']
        else:
            Results += [sciNames[item]]
c = Counter(Results).most_common()
for item in c:
    #print(item)
    if "deleted" in item[0]:
        print(item)
daa_file.close()
