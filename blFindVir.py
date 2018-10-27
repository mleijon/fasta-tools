#!/usr/bin/python3

# Something


import linecache as lica
import argparse
import sys
from fasta import FastaList


# counts the rows in infile
def countlines(infile):
    """Returns the line count."""
    with infile as f:
        for count, element in enumerate(f, 1):
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
            r = float('1'+emin)/float(e_out)
        else:
            r = float(emin)/float(e_out)
    elif float(e_in) == 0:
        e_out = 0
        r = 1
    else:
        e_out = float(e_in)
        if emin[0] == 'e':
            r = float('1' + emin)/e_out
        else:
            r = float(emin)/e_out
    return [e_out, r]


def blFindTarget():
    """Creates a dict with virus blast hits.

    Blastfile (args.b) is parsed to extract hits to viruses by
    recognizing the strings 'Virus' or 'virus' in the  'Description'
    string. Secondary hits are includeded if the ratio e-value (best hit)/
    e-value (secondary hit) is larger than args.d"""
    anchor = 'Sequences producing significant alignments: '
    blVirHit = dict()
    j = 0
    first = True
    nrOfRows = countlines(binf)
    while j <= nrOfRows:
        if anchor in lica.getline(args.b, j):
            query = lica.getline(args.b, j-6).strip()[7:]
            emin = lica.getline(args.b, j+2)[67:].split()[1].strip()
            j += 2  # Skip empty row
            hitLst = []
            if first:  # Read this only first time since constant
                offs1 = lica.getline(args.b, j).find('||') + 2
                first = False
            while lica.getline(args.b, j)[0] != '\n':
                e = lica.getline(args.b, j)[67:].split()[1].strip()
                if (args.t in lica.getline(args.b, j).casefold()
                        and ctrle(e, emin)[1] >= args.d):
                    lstEntry = dict()
                    lstEntry['Accession'] = lica.getline(args.b, j).\
                        split()[0][offs1:]
                    offs2 = offs1 + len(lstEntry['Accession'])
                    lstEntry['Description'] = lica.\
                        getline(args.b, j)[offs2:67].strip()
                    lstEntry['Bitscore'] = int(lica.getline(args.b, j)
                                               [67:].split()[0].strip())
                    lstEntry['e-value'] = ctrle(e, emin)[0]
                    lstEntry['r-e-value'] = round(ctrle(e, emin)[1], 2)
                    print(lstEntry['e-value'], lstEntry['r-e-value'], e, emin)
                    hitLst.append(lstEntry)
                j += 1
            if hitLst != []:
                blVirHit.update({query: hitLst})
        j += 1
# Returns {seqid:[{'Accession':, 'Description':,'Bitscore':,'e-value':,
# 'r-e-value':},...],...}
    return blVirHit


def crVirLst_all(blRes):
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
    filename = args.b[:args.b.find('.blast'.casefold())]+'_'+args.t\
        + '.all.txt'
    outfile = open(filename, 'w')
    for hit in sorted(VirSum, key=VirSum.__getitem__, reverse=True):
        outfile.write(hit + '|  ' + str(VirSum[hit]) + '\n')
    outfile.close()


def crVirLst_bst(blRes):
    """ Create a list of virus hits, 1 virus (the top hit)/read."""
    VirSum = dict()
    keyVS = ''
    for seqid in blRes:
        hit = blRes[seqid][0]  # Best virus hit, 1st in dict list
        keyVS = '{:12}'.format(hit['Accession']) + '|  ' + '{:55}'.format(
             hit['Description'])
        if keyVS in VirSum:
            VirSum[keyVS] += 1
        else:
            VirSum[keyVS] = 1
    hit = ''  # Hit is converted to a string
    filename = args.b[:args.b.find('.blast'.casefold())]+'_'+args.t\
        + '.bst.txt'
    outfile = open(filename, 'w')
    for hit in sorted(VirSum, key=VirSum.__getitem__, reverse=True):
        outfile.write(hit + '|  ' + str(VirSum[hit]) + '\n')
    outfile.close()


def wrTargetFa(seqidTargets):
    filename = args.b[:args.b.find('.blast'.casefold())]+'_'+args.t + '.fa'
    try:
        outfile = open(filename, 'w')
    except IOError:
        sys.exit('Output file error')
    faLst = FastaList(args.f)
    for seqid, seq in zip(faLst.id_list, faLst.seq_list):
        if seqid in seqidTargets:
            outfile.write(seq)
    outfile.close()


parser = argparse.ArgumentParser(
 description='Export blast hit containing a "target" keyword')
parser.add_argument('-f', type=str, help='fastafile', required=True)
parser.add_argument('-b', type=str, help='blastfile', required=True)
parser.add_argument('-d', type=float, default=1.0, help='sensitivity depth')
parser.add_argument('-t', type=str, default='virus', help='target')
args = parser.parse_args()
try:
    binf = open(args.b)
    finf = open(args.f)
except IOError:
    sys.exit('Input file error')
blVirHts = blFindTarget()
wrTargetFa(blVirHts.keys())
crVirLst_bst(blVirHts)
crVirLst_all(blVirHts)
