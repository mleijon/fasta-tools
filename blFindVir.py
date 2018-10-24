#!/usr/bin/python3

# Something


import linecache


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


def blFindVir(inFileName, eCut):
    anchor = 'Sequences producing significant alignments: '
    blVirHit = dict()
    j = 0
    nrOfRows = countlines(open(inFileName))
    while j <= nrOfRows:
        if anchor in linecache.getline(inFileName, j):
            query = linecache.getline(inFileName, j-6).strip()[7:]
            e = linecache.getline(inFileName, j+2)[67:].split()[1].strip()
            emin = e
            j += 2
            hitLst = []
            while linecache.getline(inFileName, j)[0] != '>':
                if (('virus' or 'Virus') in linecache.getline(inFileName, j)
                        and ctrle(e, emin)[1] > 0.1):
                    lstEntry = dict()
                    acln = len(linecache.getline(inFileName, j).split()[0]) + 1
                    lstEntry['Accession'] = linecache.getline(inFileName, j).\
                        split()[0][4:]
                    lstEntry['Description'] = linecache.\
                        getline(inFileName, j)[acln:67].strip()
                    lstEntry['Bitscore'] = int(linecache.
                                               getline(inFileName, j)[67:].
                                               split()[0].strip())

                    lstEntry['e-value'] = ctrle(e, emin)[0]
                    lstEntry['r-e-value'] = round(ctrle(e, emin)[1], 2)
                    hitLst.append(lstEntry)
                    if linecache.getline(inFileName, j+1) != '\n':
                        e = linecache.getline(inFileName, j+1)[67:].split()[1]\
                         .strip()
                j += 1
            if hitLst != []:
                blVirHit.update({query: hitLst})
        j += 1
    return blVirHit


def crVirLst_all(blRes):
    """ Create a list of virus hits. Possibly more than 1/read."""
    VirSum = dict()
    for seqid in blRes:
        for hits in blRes[seqid]:
            if hits['Accession'] + '\t' + hits['Description'] in VirSum:
                VirSum[hits['Accession'] + '\t' + hits['Description']] += 1
            else:
                VirSum[hits['Accession']+'\t' + hits['Description']] = 1
    for hit in sorted(VirSum, key=VirSum.__getitem__, reverse=True):
        print(hit, '\t', VirSum[hit],)


def crVirLst_one(blRes):
    """ Create a list of virus hits, 1 virus (top hit)/read."""
    VirSum = dict()
    hit = dict()
    for seqid in blRes:
        hit = blRes[seqid][0]
        if hit['Accession'] + '\t' + hit['Description'] in VirSum:
            VirSum[hit['Accession'] + '\t' + hit['Description']] += 1
        else:
            VirSum[hit['Accession']+'\t' + hit['Description']] = 1
    for x in sorted(VirSum, key=VirSum.__getitem__, reverse=True):
        print(x, '\t', VirSum[x],)


blVirHts = blFindVir('pool85.blast', 0.1)
crVirLst_one(blVirHts)
