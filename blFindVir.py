#!/usr/bin/python3

# Something


import linecache


# Something
def blFindVir(inFileName, eCut):
    inFile = open(inFileName)
    anchor = 'Sequences producing significant alignments: '
    hitLst = dict()
    lstEntry = dict()
    j = 0
    while linecache.getline(inFileName, j) != '':
        if linecache.getline(inFileName, j).strip() == anchor:
            eMin = float(linecache.getline(inFileName, j+1)[67:].split().
                         strip()[1])
            k = j+1
            e = float(linecache.getline(inFileName, k)[67:].split().strip()[1])
            while linecache.getline(inFileName, k)[0] != '>':
                if (('virus' or 'Virus') in linecache.getline(inFileName, k)
                        and emin/e > eCut):
                    accln = len(linecache.getline(inFileName, k).split()[0]) + 1
                    lstEntry['Accession'] = linecache.getline(inFileName, k).\
                        split()[0][5:]
                    lstEntry['Description'] = linecache.getline(inFileName, k)\
                    [accln:67].strip()

            hitLst[linecache.getline(inFileName, j-6).strip()[7:]] =\
                   linecache.getline(inFileName, k).split().strip()
            k += 1
