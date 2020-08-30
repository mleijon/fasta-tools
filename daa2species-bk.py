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
tax_classes = [ARCHAEA, BACTERIA, EUCARYOTA, VIRUSES,
               UNCLASSIFIED, OTHERSEQ]
taxsets = {CELLULAR_ORGANISM: set(), ARCHAEA: set(), BACTERIA: set(),
           EUCARYOTA: set(), VIRUSES: set(), UNCLASSIFIED: set(),
           OTHERSEQ: set()}

# THese files are needed and downloaded from ncbi:
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

names_file = "names.dmp"
# Taxonomy names file has these fields:
#
# tax_id					-- the id of node associated with this name
# name_txt		     		-- name itself
# unique name				-- the unique variant of this name if name not
#                              unique
# name class				-- (synonym, common name, ...)

deleted_file = 'delnodes.dmp'
# Deleted nodes (nodes that existed but were deleted) file field:
#
# tax_id					-- deleted node id

merged_file = "merged.dmp"
# Merged nodes file fields:
#
# old_tax_id                -- id of nodes which has been merged
# new_tax_id                -- id of nodes which is result of merging


lineage_file = "taxidlineage.dmp"
# Taxonomy id lineage file fields:
#
# tax_id                    -- node id
# lineage                   -- sequence of node ids separated by space denoting
#                              nodes' ancestors starting from the most distant
#                              one and ending with the immediate one


# Count the sequence reads (rows) in the input diamond file (.daa)
def count_rows(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# Counts the number of reads without hits, '0' in .daa file, and the number of
# reads assigned to 'ROOT', '1' in .daa file. The interpretation of '1' as
# a taxonomy id is that the reads have hit to the database but of equal
# similarity to enrries in the database of widely different taxa.
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


# Create sets with the taxonomic classes: ARCHAEA, BACTERIA, EUCARYOTA, VIRUSES,
# UNCLASSIFIED, OTHERSEQ
def create_taxsets(fname):
    with open(fname) as f:
        for line in f:
            line = line.replace('|', ' ')
            item = list(map(str.strip, line.split()[0:3]))
            for taxclass in tax_classes:
                if taxclass in item:
                    taxsets[taxclass].add(item[0])
                    break
    return


def count_txclass(fname, taxclass):
    if taxclass == ARCHAEA:
        taxclass_set = taxsets[ARCHAEA]
    elif taxclass == BACTERIA:
        taxclass_set = taxsets[BACTERIA]
    elif taxclass == EUCARYOTA:
        taxclass_set = taxsets[EUCARYOTA]
    elif taxclass == VIRUSES:
        taxclass_set = taxsets[VIRUSES]
    elif taxclass == UNCLASSIFIED:
        taxclass_set = taxsets[UNCLASSIFIED]
    elif taxclass == OTHERSEQ:
        taxclass_set = taxsets[OTHERSEQ]
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


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input diamond file',
                        required=True)
    PARSER.add_argument('-s', action='store_true',
                        help='switch for summary file output')
    ARGS = PARSER.parse_args()

    sciNames = dict()
    mergedTaxids = dict()
    deletedTaxids = set()
    Results = []

    nr_of_reads = count_rows(ARGS.f)
    nr_of_nohits = count_0and1(ARGS.f)
    create_taxsets("taxidlineage.dmp")

    # Creates dictionary with merged taxids
    with open(merged_file) as f:
        for line in f:
            item = line.split('|')
            mergedTaxids.update({item[0].strip(): item[1].strip()})

    # Creates a set with deleted taxids
    with open(deleted_file) as f:
        for line in f:
            item = line.split('|')[0].strip()
            deletedTaxids.add(item)

    # Creates a dictionary with key:values equal taxid:scientific name
    with open(names_file) as f:
        for line in f:
            item = line.split('|')
            if "scientific" in item[3]:
                sciNames.update({item[0].strip(): item[1].strip()})

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
    print('\nNo hits to db: {:>35} ({} %)'.format(nr_of_nohits[0], (
        round(100 * (nr_of_nohits[0] / nr_of_reads), 2))))
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
    print('UNCLASSIFIED: {:>36} ({} %)'.format(
        count_txclass(ARGS.f, UNCLASSIFIED), (round(100 * (
                count_txclass(ARGS.f, UNCLASSIFIED) / nr_of_reads), 2))))
    print('OTHERSEQ: {:>40} ({} %)'.format(count_txclass(ARGS.f, OTHERSEQ), (
        round(100 * (count_txclass(ARGS.f, OTHERSEQ) / nr_of_reads), 2))))
