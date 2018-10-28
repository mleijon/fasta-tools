#!/usr/bin/python3
"""Something"""


import linecache as lica
import argparse
import sys
from fasta import FastaList


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
    """Catches blast file eceptions.

    The blast file of Decypher loses the digit (x) before the 'e' in the
    e-value (i.e.) xe-yyy when yyy is larger than 099. '1' is added to fix
    this. Division by zero, if the e-value is zero is also avoided by this
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


def blFindTarget():
    """Creates a dict with virus blast hits.

    Blastfile (ARGS.b) is parsed to extract hits to viruses by
    recognizing the strings 'Virus' or 'virus' in the  'Description'
    string. Secondary hits are includeded if the ratio e-value (best hit)/
    e-value (secondary hit) is larger than ARGS.d"""
    anchor = 'Sequences producing significant alignments: '
    blVirHit = dict()
    j = 0
    first = True
    nrOfRows = count_li(BL_IN)
    while j <= nrOfRows:
        if anchor in lica.getline(ARGS.b, j):
            query = lica.getline(ARGS.b, j-6).strip()[7:]
            emin = lica.getline(ARGS.b, j+2)[67:].split()[1].strip()
            j += 2  # Skip empty row
            hitLst = []
            if first:  # Read this only first time since constant
                offs1 = lica.getline(ARGS.b, j).find('||') + 2
                first = False
            while lica.getline(ARGS.b, j)[0] != '\n':
                e_val = lica.getline(ARGS.b, j)[67:].split()[1].strip()
                if (ARGS.t in lica.getline(ARGS.b, j).casefold()
                        and ctrle(e_val, emin)[1] >= ARGS.d):
                    lstEntry = dict()
                    lstEntry['Accession'] = lica.getline(ARGS.b, j).\
                        split()[0][offs1:]
                    offs2 = offs1 + len(lstEntry['Accession'])
                    lstEntry['Description'] = lica.\
                        getline(ARGS.b, j)[offs2:67].strip()
                    lstEntry['Bitscore'] = int(lica.getline(ARGS.b, j)
                                               [67:].split()[0].strip())
                    lstEntry['e-value'] = ctrle(e_val, emin)[0]
                    lstEntry['r-e-value'] = round(ctrle(e_val, emin)[1], 2)
                    hitLst.append(lstEntry)
                j += 1
            if hitLst != []:
                blVirHit.update({query: hitLst})
        j += 1
# Returns {seqid:[{'Accession':, 'Description':,'Bitscore':,'e-value':,
# 'r-e-value':},...],...}
    return blVirHit


def wr_deep_tar(blRes):
    """ Create a list of virus hits. Possibly more than 1/read."""
    tar_sum = dict()
    tar_sum_key = ''
    for seqid in blRes:
        for hit in blRes[seqid]:
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
    outfile.write('{:15}'.format('Accession') + '{:60}'.format('Description')
                  + '{:10}'.format('e-min') + '{:5}'.format('Nr-of-reads')
                  + '\n\n')
    for hits in sorted(tar_sum, key=lambda hits: tar_sum[hits]['readSum'],
                       reverse=True):
        outfile.write('{:12}'.format(hits) + '|  ' + '{:55}'.format(
            tar_sum[hits]['Description']) + '|  ' + '{:8}'.format(
                tar_sum[hits]['emin']) + '  |  ' + '{:7}'.format(
                    tar_sum[hits]['readSum']) + '\n')
    outfile.close()


def wr_top_tar(blRes):
    """ Create a list of target hits.

    If depth is set to 1 (default) Only hits where the target are the best
    are give (i.e. no other organism have a better hit). If the depth i set
    lower than 1, the best virus hit within the limit of the depth threshold
    is given.
    """
    tar_sum = dict()
    for seqid in blRes:
        hit = blRes[seqid][0]  # Best target hit, 1st in dict list
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
    outfile.write('{:15}'.format('Accession') + '{:60}'.format('Description')
                  + '{:10}'.format('e-min') + '{:5}'.format('Nr-of-reads')
                  + '\n\n')
    for hits in sorted(tar_sum, key=lambda hits: tar_sum[hits]['readSum'],
                       reverse=True):
        outfile.write('{:12}'.format(hits) + '|  ' + '{:55}'.format(
            tar_sum[hits]['Description']) + '|  ' + '{:8}'.format(
                tar_sum[hits]['emin']) + '  |  ' + '{:7}'.format(
                    tar_sum[hits]['readSum']) + '\n')
    outfile.close()


def wr_fa_tar(seqidTargets):
    """Writes a fasta file with the target hitting sequences"""
    filename = ARGS.b[:ARGS.b.find('.blast'.casefold())]+'_'+ARGS.t + '.fa'
    try:
        outfile = open(filename, 'w')
    except IOError:
        sys.exit('Output file error')
    faLst = FastaList(ARGS.f)
    for seqid, seq in zip(faLst.id_list, faLst.seq_list):
        if seqid in seqidTargets:
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
TAR_HITS = blFindTarget()
wr_fa_tar(TAR_HITS.keys())
wr_top_tar(TAR_HITS)
wr_deep_tar(TAR_HITS)
print('Done!      ')
