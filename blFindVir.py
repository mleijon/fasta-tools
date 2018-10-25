#!/usr/bin/python3

# Something


import linecache
import argparse
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
            r = float('1' + emin)/float(e_in)
        else:
            r = float(emin)/float(e_in)
    return [e_out, r]


def blFindTarget(arg.b, arg.d, arg.t):
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
        if anchor in linecache.getline(args.b, j):
            query = linecache.getline(args.b, j-6).strip()[7:]
            e = linecache.getline(args.b, j+2)[67:].split()[1].strip()
            emin = e  # First hit = best hit
            j += 2  # Skip empty row
            hitLst = []
            if first:  # Read this only first time since constant
                offs1 = linecache.getline(args.b, j).find('||') + 2
                first = False
            while linecache.getline(args.b, j)[0] != '\n':
                if (args.t in linecache.getline(args.b, j).casefold()
                        and ctrle(e, emin)[1] > args.d):
                    lstEntry = dict()
                    lstEntry['Accession'] = linecache.getline(args.b, j).\
                        split()[0][offs1:]
                    offs2 = offs1 + len(lstEntry['Accession'])
                    lstEntry['Description'] = linecache.\
                        getline(args.b, j)[offs2:67].strip()
                    lstEntry['Bitscore'] = int(linecache.getline(args.b, j)
                                               [67:].split()[0].strip())
                    lstEntry['e-value'] = ctrle(e, emin)[0]
                    lstEntry['r-e-value'] = round(ctrle(e, emin)[1], 2)
                    hitLst.append(lstEntry)
                    if linecache.getline(args.b, j+1) != '\n':
                        e = linecache.getline(args.b, j+1)[67:].split()[1]\
                         .strip()
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
    for seqid in blRes:
        for hit in blRes[seqid]:
            if hit['Accession'] + '\t' + hit['Description'] in VirSum:
                VirSum[hit['Accession'] + '\t' + hit['Description']] += 1
            else:
                VirSum[hit['Accession']+'\t' + hit['Description']] = 1
    hit = ''  # Here hit is converted to a string
    for hit in sorted(VirSum, key=VirSum.__getitem__, reverse=True):
        print(hit, '\t', VirSum[hit],)


def crVirLst_bst(blRes):
    """ Create a list of virus hits, 1 virus (the top hit)/read."""
    VirSum = dict()
    for seqid in blRes:
        hit = blRes[seqid][0]  # Best virus hit, 1st in dict list
        if hit['Accession'] + '\t' + hit['Description'] in VirSum:
            VirSum[hit['Accession'] + '\t' + hit['Description']] += 1
        else:
            VirSum[hit['Accession']+'\t' + hit['Description']] = 1
    hit = ''  # Here hit is converted to a string
    for hit in sorted(VirSum, key=VirSum.__getitem__, reverse=True):
        print(hit, '\t', VirSum[hit],)

def wrTargetFa():
    pass

parser = argparse.ArgumentParser(
description='Export blast hit containing a "target" keyword')
parser.add_argument('-f', type=str, help='fastafile')
parser.add_argument('-b', type=str, help='blastfile')
parser.add_argument('-d', type=float, default=0.1, help='sensitivity depth')
parser-add_argument('t', type=str, default='virus',help='target')
args = parser.parse_args()
try:
    binf = open(args.b)
    finf = open(args.f)
except IOError:
    sys.exit('Input file error')
blVirHts = blFindTarget(arg.b, arg.d, arg.t)
seqidTgtLst = blFindTarget.keys()
faLst = FastaList(args.f)
crVirLst_bst(blVirHts)
crVirLst_all(blVirHts)
