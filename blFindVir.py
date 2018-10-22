#!/usr/bin/python3

# Something


import linecache


# counts the rows in infile
def countlines(infile):
    with infile as f:
        for count, element in enumerate(f, 1):
            pass
    infile.close()
    return count


# Something
def ctrle(e_in, emin):
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
    hitLst = []
    lstEntry = dict()
    blVirHit = dict()
    j = 0
    nrOfRows = countlines(open(inFileName))
    while j <= nrOfRows:
        if anchor in linecache.getline(inFileName, j):
            query = linecache.getline(inFileName, j-6).strip()[7:]
            e = linecache.getline(inFileName, j+2)[67:].split()[1].strip()
            emin = e
            j += 2
            while linecache.getline(inFileName, j)[0] != '>':
                if (('virus' or 'Virus') in linecache.getline(inFileName, j)
                        and ctrle(e, emin)[1] > 0.1):
                    acln = len(linecache.getline(inFileName, j).split()[0]) + 1
                    lstEntry['Accession'] = linecache.getline(inFileName, j).\
                        split()[0][5:]
                    lstEntry['Description'] = linecache.\
                        getline(inFileName, j)[acln:67].strip()
                    lstEntry['Bitscore'] = int(linecache.
                                               getline(inFileName, j)[67:].
                                               split()[0].strip())

                    lstEntry['e-value'] = ctrle(e, emin)[0]
                    lstEntry['r-e-value'] = round(ctrle(e, emin)[1], 2)
                    print(lstEntry)
                    hitLst.append(lstEntry)
                    if linecache.getline(inFileName, j+1) != '\n':
                        e = linecache.getline(inFileName, j+1)[67:].split()[1]\
                         .strip()
                j += 1
            blVirHit[query] = hitLst
        j += 1

    return blVirHit


blFindVir('2C.fastq.blast', 0.1)
