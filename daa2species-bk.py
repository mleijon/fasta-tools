#!/usr/bin/python3

from collections import Counter, defaultdict
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
taxcounts = defaultdict(int)

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


def count_txclass(fname):
    with open(fname) as f:
        for row in f:
            taxid = row.split()[1].strip()
            if taxid in taxsets[ARCHAEA]:
                taxcounts[ARCHAEA] += 1
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[ARCHAEA]):
                taxcounts[ARCHAEA] += 1
            elif taxid in taxsets[BACTERIA]:
                taxcounts[BACTERIA] += 1
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[BACTERIA]):
                taxcounts[BACTERIA] += 1
            elif taxid in taxsets[EUCARYOTA]:
                taxcounts[EUCARYOTA] += 1
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[EUCARYOTA]):
                taxcounts[EUCARYOTA] += 1
            elif taxid in taxsets[VIRUSES]:
                taxcounts[VIRUSES] += 1
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[VIRUSES]):
                taxcounts[VIRUSES] += 1
            elif taxid in taxsets[UNCLASSIFIED]:
                taxcounts[UNCLASSIFIED] += 1
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[UNCLASSIFIED]):
                taxcounts[UNCLASSIFIED] += 1
            elif taxid in taxsets[OTHERSEQ]:
                taxcounts[OTHERSEQ] += 1
            elif (taxid in mergedTaxids.keys()) and (mergedTaxids[taxid]
                                                     in taxsets[OTHERSEQ]):
                taxcounts[OTHERSEQ] += 1
    return


def parse_files(m, d, n):
    # Creates dictionary with merged taxids
    with open(m) as f:
        for line in f:
            item = line.split('|')
            mergedTaxids.update({item[0].strip(): item[1].strip()})

    # Creates a set with deleted taxids
    with open(d) as f:
        for line in f:
            item = line.split('|')[0].strip()
            deletedTaxids.add(item)

    # Creates a dictionary with key:values equal taxid:scientific name
    with open(n) as f:
        for line in f:
            item = line.split('|')
            if "scientific" in item[3]:
                sciNames.update({item[0].strip(): item[1].strip()})


def write_summary (fname):
    with open(fname, 'w') as f:
        f.write('Total nr of reads: {:>31}\n'.format(nr_of_reads))
        f.write('\nNo hits to db: {:>35} ({} %)\n'.format(nr_of_nohits[0], (
            round(100 * (nr_of_nohits[0] / nr_of_reads), 2))))
        f.write('Could not be classified by db: {:>19} ({} %)\n'.format(
            nr_of_nohits[1], (round(100 * (nr_of_nohits[1] / nr_of_reads), 2))))
        f.write('Taxid deleted from db: {:>27} ({} %)\n'.format(
            deleted_count, (round(100 * (deleted_count / nr_of_reads), 2))))
        f.write('\nClassified by db: {:>32} ({} %)\n'.format(
            SUM_IN_DB, (round(100 * (SUM_IN_DB / nr_of_reads), 2))))
        f.write('------------------------------------------------------------\n')
        f.write('CELLULAR ORGANISM: {:>31} ({} %)\n'.format(cellorg_count, (
            round(100 * (cellorg_count / nr_of_reads), 2))))
        f.write('ARCHAEA: {:>41} ({} %)\n'.format(taxcounts[ARCHAEA], (
            round(100 * (taxcounts[ARCHAEA] / nr_of_reads), 2))))
        f.write('BACTERIA: {:>40} ({} %)\n'.format(taxcounts[BACTERIA], (
            round(100 * (taxcounts[BACTERIA] / nr_of_reads), 2))))
        f.write('EUCARYOTA: {:>39} ({} %)\n'.format(taxcounts[EUCARYOTA], (
            round(100 * (taxcounts[EUCARYOTA] / nr_of_reads), 2))))
        f.write('VIRUS: {:>43} ({} %)\n'.format(taxcounts[VIRUSES], (
            round(100 * (taxcounts[VIRUSES] / nr_of_reads), 2))))
        f.write('UNCLASSIFIED: {:>36} ({} %)\n'.format(
            taxcounts[UNCLASSIFIED], (round(100 * (
                    taxcounts[UNCLASSIFIED] / nr_of_reads), 2))))
        f.write('OTHERSEQ: {:>40} ({} %)\n'.format(taxcounts[OTHERSEQ], (
            round(100 * (taxcounts[OTHERSEQ] / nr_of_reads), 2))))


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
    parse_files(merged_file, deleted_file, names_file)
    count_txclass(ARGS.f)

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

    SUM_IN_DB = sum(taxcounts.values()) + cellorg_count
    if ARGS.s:
        summary_file_name = ARGS.f.split('.')[0] + '_summary.txt'
        write_summary(summary_file_name)
