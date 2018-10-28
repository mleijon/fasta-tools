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
    VirSum = dict()
    keyVS = ''
    for seqid in blRes:
        for hit in blRes[seqid]:
            keyVS = '{:12}'.format(hit['Accession']) + '|  ' + '{:55}'.format(
                hit['Description'])
            if keyVS in VirSum:
                VirSum[keyVS] += 1
            else:
                VirSum[keyVS] = 1
    hit = ''  # Hit is converted to a string
    filename = ARGS.b[:ARGS.b.find('.blast'.casefold())]+'_'+ARGS.t\
        + '.deep.txt'
    outfile = open(filename, 'w')
    for hit in sorted(VirSum, key=VirSum.__getitem__, reverse=True):
        outfile.write(hit + '|  ' + str(VirSum[hit]) + '\n')
    outfile.close()


def wr_top_tar(blRes):
    """ Create a list of virus hits, 1 virus (the top hit)/read."""
    VirSum = dict()
    for seqid in blRes:
        hit = blRes[seqid][0]  # Best target hit, 1st in dict list
        keyVS = hit['Accession']
        VirSum.setdefault(keyVS, {'Description': hit['Description'],
                                  'emin': 1, 'readSum': 0})
        if keyVS in VirSum.keys():
            VirSum[keyVS]['readSum'] += 1
            if hit['e-value'] < VirSum[keyVS]['emin']:
                VirSum[keyVS]['emin'] = hit['e-value']
    filename = ARGS.b[:ARGS.b.find('.blast'.casefold())]+'_'+ARGS.t\
        + '.top.txt'
    outfile = open(filename, 'w')
    for hit in sorted(VirSum, key=lambda hit: VirSum[hit]['readSum'],
                      reverse=True):
        outfile.write(hit + '|  ' + str(VirSum[hit]) + '\n')
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
