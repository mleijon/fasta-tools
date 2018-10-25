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


def blFindTarget(inFileName, eCut, target):
    """Creates a dict with virus blast hits.

    Blastfile (inFileName) is parsed to extract hits to viruses by
    recognizing the strings 'Virus' or 'virus' in the  'Description'
    string. Secondary hits are includeded if the ratio e-value (best hit)/
    e-value (secondary hit) is larger than eCut"""
    anchor = 'Sequences producing significant alignments: '
    blVirHit = dict()
    j = 0
    first = True
    nrOfRows = countlines(open(inFileName))
    while j <= nrOfRows:
        if anchor in linecache.getline(inFileName, j):
            query = linecache.getline(inFileName, j-6).strip()[7:]
            e = linecache.getline(inFileName, j+2)[67:].split()[1].strip()
            emin = e  # First hit = best hit
            j += 2  # Skip empty row
            hitLst = []
            if first:  # Read this only first time since constant
                offs1 = linecache.getline(inFileName, j).find('||') + 2
                first = False
            while linecache.getline(inFileName, j)[0] != '\n':
                if (virus in linecache.getline(inFileName, j).casefold()
                        and ctrle(e, emin)[1] > eCut):
                    lstEntry = dict()
                    lstEntry['Accession'] = linecache.getline(inFileName, j).\
                        split()[0][offs1:]
                    offs2 = offs1 + len(lstEntry['Accession'])
                    lstEntry['Description'] = linecache.\
                        getline(inFileName, j)[offs2:67].strip()
                    lstEntry['Bitscore'] = int(linecache.getline(inFileName, j)
                                               [67:].split()[0].strip())
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
    return blVirHit # {seqid:[{'Accession':, 'Description':,'Bitscore':,
                    # 'e-value':, 'r-e-value':},...],...}


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


blVirHts = blFindTarget('pool85.blast', 0.1,'virus')
crVirLst_bst(blVirHts)
