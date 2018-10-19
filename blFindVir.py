#!/usr/bin/python3

# Something


import linecache


# Something
def blFindVir(inFileName, eCut):
    anchor = 'Sequences producing significant alignments: '
    hitLst = []
    lstEntry = dict()
    blVirHit = dict()
    j = 0
    while linecache.getline(inFileName, j) != ' ':
        if anchor in linecache.getline(inFileName, j):
            print(linecache.getline(inFileName, j+1)[67:])
            eMin = float(linecache.getline(inFileName, j+2)[67:].split()[1].
                         strip())
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
                    lstEntry['e-value'] = float(linecache.
                                                getline(inFileName, j)[67:].
                                                split()[1].strip())
                    lstEntry['r-e-value'] = eMin/e
                    hitLst.append(lstEntry)
                    if linecache.getline(inFileName, j+1)[0] != '>':
                        print(linecache.getline(inFileName, j+1))
                        e = float(linecache.getline(inFileName, j+1)[67:].
                                  split()[1].strip())
                    j += 1
            blVirHit[linecache.getline(inFileName, j-6).strip()[7:]] = hitLst
        j += 1


blFindVir('2C.fastq.blast', 0.1)
