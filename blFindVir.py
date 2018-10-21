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
            if linecache.getline(inFileName, j+2)[67:].split()[1].strip()[0]\
               == 'e':
                    eMin = float('1' + linecache.getline(inFileName, j+2)[67:].
                                 split()[1].strip())
            else:
                eMin = float(linecache.getline(inFileName, j+2)[67:].
                             split()[1].strip())
            j += 2
            e = eMin
            while linecache.getline(inFileName, j)[0] != '>':
                if (('virus' or 'Virus') in linecache.getline(inFileName, j)
                        and eMin/e > eCut):
                    acln = len(linecache.getline(inFileName, j).split()[0]) + 1
                    lstEntry['Accession'] = linecache.getline(inFileName, j).\
                        split()[0][5:]
                    lstEntry['Description'] = linecache.\
                        getline(inFileName, j)[acln:67].strip()
                    lstEntry['Bitscore'] = int(linecache.
                                               getline(inFileName, j)[67:].
                                               split()[0].strip())
                    if linecache.getline(inFileName, j)[67:].split()[1].\
                       strip()[0] == 'e':
                        lstEntry['e-value'] = float('1' + linecache.
                                                    getline(inFileName, j)
                                                    [67:].split()[1].strip())
                    else:
                        lstEntry['e-value'] = float(linecache.
                                                    getline(inFileName, j)
                                                    [67:].split()[1].strip())
                    lstEntry['r-e-value'] = eMin/e
                    hitLst.append(lstEntry)
                    if linecache.getline(inFileName, j+1)[0] != '\n':
                        if linecache.getline(inFileName, j+1)[67:].split()[1].\
                         strip()[0] == 'e':
                            e = float('1' + linecache.getline(inFileName, j+1)
                                      [67:].split()[1].strip())
                        else:
                            e = float(linecache.getline(inFileName, j+1)
                                      [67:].split()[1].strip())
                j += 1
            blVirHit[query] = hitLst
        j += 1

    return blVirHit


print(blFindVir('newf.blast', 0.1)["NODE_26_length_7842_cov_3.707842"][1]["Accession"])
