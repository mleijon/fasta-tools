#!/usr/bin/python3
"""Finds the blast hits containing a target string in the description,
irrespective case.

Creates three output files.
"""


import linecache as lica
import argparse
import sys


# counts the rows in infile
def count_li(infile):
    """Returns the line count."""
    count = 1
    with infile as file:
        for count, element in enumerate(file, start=1):
            pass
    infile.close()
    return count


# Catches blast file exceptions
def ctrle(e_in, emin):
    """Catches blast file exceptions and cakculated e-value ratios.

    1) The blast file of Decypher loses the digit (x) before the 'e' in the
    e-value (i.e.) xe-yyy when yyy is larger than 099. '1' is added to fix
    this. 2) Division by zero, if the e-value is zero is also handeled by this
    function."""
    if e_in[0] == 'e':
        e_out = float('1' + e_in)
        if emin[0] == 'e':
            ratio = float('1'+emin)/float(e_out)
        else:
            ratio = float(emin)/float(e_out)
    elif float(e_in) == 0:
        e_out = 0
        ratio = 1
    else:
        e_out = float(e_in)
        if emin[0] == 'e':
            ratio = float('1' + emin)/e_out
        else:
            ratio = float(emin)/e_out
    return [e_out, ratio]


def find_targets():
    """Creates a dict with virus blast hits.

    Blastfile (ARGS.b) is parsed to extract hits to a target by
    recognizing the string 'target', independent of case, in the 'Description'
    string. Secondary hits are includeded if the ratio e-value (best hit)/
    e-value (secondary hit) is larger than ARGS.d
    """
    anchor = 'Sequences producing significant alignments: '
    target_hits = dict()
    j = 0
    first = True
    nr_rows = count_li(BL_IN)
    while j <= nr_rows:
        if anchor in lica.getline(ARGS.b, j):
            query = lica.getline(ARGS.b, j-6).strip()[7:]
            emin = lica.getline(ARGS.b, j+2)[67:].split()[1].strip()
            j += 2  # Skip empty row
            tar_lst = []
            if first:  # Read this only first time since constant
                offs1 = lica.getline(ARGS.b, j).find('||') + 2
                first = False
            while lica.getline(ARGS.b, j)[0] != '\n':
                e_val = lica.getline(ARGS.b, j)[67:].split()[1].strip()
                if (ARGS.t.casefold() in lica.getline(ARGS.b, j).casefold()
                        and ctrle(e_val, emin)[1] >= ARGS.d):
                    tar_hit = dict()
                    tar_hit['Accession'] = lica.getline(ARGS.b, j).\
                        split()[0][offs1:]
                    offs2 = offs1 + len(tar_hit['Accession'])
                    tar_hit['Description'] = lica.\
                        getline(ARGS.b, j)[offs2:67].strip()
                    tar_hit['Bitscore'] = int(lica.getline(ARGS.b, j)
                                              [67:].split()[0].strip())
                    tar_hit['e-value'] = ctrle(e_val, emin)[0]
                    tar_hit['r-e-value'] = round(ctrle(e_val, emin)[1], 2)
                    tar_lst.append(tar_hit)
                j += 1
            if tar_lst:
                target_hits.update({query: tar_lst})
        j += 1
# Returns {seqid:[{'Accession':, 'Description':,'Bitscore':,'e-value':,
# 'r-e-value':},...],...}
    return target_hits


def wr_deep_tar(bl_result):
    """ Create a list of target hits.

    If depth is set to 1 (default) All hits where the targets are the best
    are given (i.e. no other organism have a better hit). If the depth is set
    lower than 1, all target hits within the limit of the depth threshold
    are given.
    """
    tar_sum = dict()
    for seqid in bl_result:
        for hit in bl_result[seqid]:
            tar_sum_key = hit['Accession']
            tar_sum.setdefault(tar_sum_key, {'Description': hit['Description'],
                                             'emin': 1, 'readSum': 0})
            if tar_sum_key in tar_sum.keys():
                tar_sum[tar_sum_key]['readSum'] += 1
                if hit['e-value'] < tar_sum[tar_sum_key]['emin']:
                    tar_sum[tar_sum_key]['emin'] = hit['e-value']
    filename = ARGS.b[:ARGS.b.find('.blast'.casefold())]+'_'+ARGS.t\
        + '.deep.txt'
    outfile = open(filename, 'w')
    outfile.write('{:19}'.format('Accession') + '{:60}'.format('Description')
                  + '{:10}'.format('e-min') + '{:5}'.format('Nr-of-reads')
                  + '\n\n')
    for hits in sorted(tar_sum, key=lambda hits: tar_sum[hits]['readSum'],
                       reverse=True):
        outfile.write('{:16}'.format(hits) + '|  ' + '{:55}'.format(
            tar_sum[hits]['Description']) + '|  ' + '{:8}'.format(
                tar_sum[hits]['emin']) + '  |  ' + '{:7}'.format(
                    tar_sum[hits]['readSum']) + '\n')
    outfile.close()


def wr_top_tar(bl_result):
    """ Create a list of target hits.

    If depth is set to 1 (default) Only hits where hits to the target are the
    best are given (i.e. no other organism have a better hit). If the depth is
    set lower than 1, the best target hit within the limit of the depth
    threshold is given.
    """
    tar_sum = dict()
    for seqid in bl_result:
        hit = bl_result[seqid][0]  # Best target hit, 1st in dict list
        tar_sum_key = hit['Accession']
        tar_sum.setdefault(tar_sum_key, {'Description': hit['Description'],
                                         'emin': 1, 'readSum': 0})
        if tar_sum_key in tar_sum.keys():
            tar_sum[tar_sum_key]['readSum'] += 1
            if hit['e-value'] < tar_sum[tar_sum_key]['emin']:
                tar_sum[tar_sum_key]['emin'] = hit['e-value']
    filename = ARGS.b[:ARGS.b.find('.blast'.casefold())]+'_'+ARGS.t\
        + '.top.txt'
    outfile = open(filename, 'w')
    outfile.write('{:19}'.format('Accession') + '{:60}'.format('Description')
                  + '{:10}'.format('e-min') + '{:5}'.format('Nr-of-reads')
                  + '\n\n')
    for hits in sorted(tar_sum, key=lambda hits: tar_sum[hits]['readSum'],
                       reverse=True):
        outfile.write('{:16}'.format(hits) + '|  ' + '{:55}'.format(
            tar_sum[hits]['Description']) + '|  ' + '{:8}'.format(
                tar_sum[hits]['emin']) + '  |  ' + '{:7}'.format(
                    tar_sum[hits]['readSum']) + '\n')
    outfile.close()


def wr_fa_tar(seq_ids):
    """Writes fasta-file with sequences that matches the target"""

    from fasta import FastaList
    filename = ARGS.b[:ARGS.b.find('.blast'.casefold())]+'_'+ARGS.t + '.fa'
    try:
        outfile = open(filename, 'w')
    except IOError:
        sys.exit('Output file error')
    fa_lst = FastaList(ARGS.f)
    for seqid, seq in zip(fa_lst.id_list, fa_lst.seq_list):
        if seqid in seq_ids:
            outfile.write(seq)
    outfile.close()


PARSER = argparse.ArgumentParser(description='Export blast hit containing a\
                                     target keyword')
PARSER.add_argument('-f', type=str, help='fastafile', required=True)
PARSER.add_argument('-b', type=str, help='blastfile', required=True)
PARSER.add_argument('-d', type=float, default=1.0, help='sensitivity depth')
PARSER.add_argument('-t', type=str, default='virus', help='target')
ARGS = PARSER.parse_args()
try:
    BL_IN = open(ARGS.b)
    FA_IN = open(ARGS.f)
except IOError:
    sys.exit('Input file error')
print('Parsing...', end='\r', flush=True)
TAR_HITS = find_targets()
wr_fa_tar(TAR_HITS.keys())
wr_top_tar(TAR_HITS)
wr_deep_tar(TAR_HITS)
print('Done!      ')
