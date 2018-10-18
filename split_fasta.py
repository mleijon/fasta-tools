#!/usr/bin/python3

# Creates a list of dictionaries of the sequence parameters assuming SPADES
# (http://cab.spbu.ru/software/spades/) contigs. NODE (the sequence number),
# length (length of sequence), cov (read coverage of the contig).
import linecache
import argparse
import sys
import os
from fasta import FastaList


class spadesFa(FastaList):
    """Class derived from SPades contig fasta file.

       The SPades fasta-file header lines contain information on read coverage
       and contig length and a generic "NODE" number. These paremeters are
       stored in a list of dictionaries as a property of the class. Other
       proprties are nr of contigs, the average coverage, the average contig
       length and a boolean parameter being true if the fasta file is generated
       by SPades (http://cab.spbu.ru/software/spades/). The class has a
       function writing the parameters to a file. The class inherits the
       FastaList class and uses the id_list propery of this class in the
       WrFiles function of the main program."""

    def __init__(self, spades_file):
        super(spadesFa, self).__init__(args.f)
        sumCov = 0
        sumLen = 0
        t = []
        seqPar = dict()
        self.seqParLst = []
        self.is_spadesFa = True
        for line in spades_file:
            if line.startswith('>'):
                line = line.strip()
                t = line.split('_')
                if not (t[0][1:] == 'NODE' and t[2] == 'length'
                        and t[4] == 'cov'):
                    self.is_spadesFa = False
                    break
                seqPar = {t[0][1:]: int(t[1]), t[2]: int(t[3]), t[4]:
                          float(t[5])}  # parses the id strings of fa-file
                self.seqParLst.append(seqPar)
        if self.is_spadesFa:
            self.nrOfContig = len(self.seqParLst)
            for i in range(len(self.seqParLst)):
                sumCov += self.seqParLst[i]['cov']
                sumLen += self.seqParLst[i]['length']
            self.avLen = round(sumLen/self.nrOfContig, 0)
            self.avCov = round(sumCov/self.nrOfContig, 0)

# Writes a parameter file based on the ID-string created by Spades and some
# other parameters
    def wrPar2File(self, parfile, inpname):
            parfile.write('Parameters for: '+inpname+'\n\n')
            parfile.write('Number of contigs: %d\n' % self.nrOfContig)
            parfile.write('Average contig length: %d\n' % self.avLen)
            parfile.write('Average contig coverage: %d\n' % self.avCov)
            parfile.write('_________________________________________\n\n')


class blastTbl(object):
    """Create a class from a blast table from the corrsponding contigs created
       by SPades.

       Properties are a lst containing fasta ids for contigs with no blast hits
       and the number of sequences without blast hits"""

    def __init__(self, inFile):
        super(blastTbl, self).__init__()
        self.noHits = self.noBlastHit(inFile)
        self.nrNoHits = len(self.noHits)

# Creates dictionary of blast *no hits* ids
    def noBlastHit(self, inFile):
        no_hit = '***** No hits found ******'
        no_hits = dict()
        j = 0
        for line in inFile:
            j += 1
            if line.strip() == no_hit:
                no_hits[linecache.getline(args.b, j-6).strip()[7:]] = j-6
        linecache.clearcache()
        return no_hits

# Write filtered fasta file not including sequeces with no blast hits and a
# fasta file including exlusively these contigs. A parameter file is also
# written using the FastaList class


def WrFiles():
    if args.f[-2:] == ('gz' or 'GZ'):
        fahits = 'hits_'+args.f[:-3]
        fanohits = 'nohits_'+args.f[:-3]
    else:
        fahits = 'hits_'+args.f
        fanohits = 'nohits_'+args.f
    if fahits[-2:] == ('fq' or 'FQ'):
        fahits = fahits[:-3]
        fanohits = fanohits[:-3]
    elif fahits[-5:] == ('fastq' or 'FASTQ'):
        fahits = fahits[:-5]
        fanohits = fanohits[:-5]
    try:
        FiHit = open(fahits+'fa', 'w')
        FiNoHit = open(fanohits+'fa', 'w')
    except IOError:
        sys.exit('error writing file')
    j = 0
    for ids in spFaLst.id_list:
        if ids not in blLst.noHits:
            FiHit.write(spFaLst.seq_list[j])
        else:
            FiNoHit.write(spFaLst.seq_list[j])
        j += 1
    FiHit.close()
    FiNoHit.close()
    if args.p and spFaLst.is_spadesFa:
        if '.' in args.f:
            name = args.f[:args.f.find('.')]+'.par'
            fipa = open(name, 'w')
        else:
            fipa = open(args.f+'.par', 'w')
        spFaLst.wrPar2File(fipa, args.f)
        fihi = open(fahits)
        faHitLst = spadesFa(fihi)
        faHitLst.wrPar2File(fipa, fahits)
        percent_remove = round(100*blLst.nrNoHits/len(spFaLst.seq_list), 0)
        fipa.write('%d sequences removed (%g%%)\n' % (blLst.nrNoHits,
                                                      percent_remove))
        fipa.close()
        fihi.close()
    elif args.p:
        fipa = open('par_'+args.f[:args.f.find('.')], 'w')
        fipa.write('No parameters, not contigs fasta file created by spades.')
        fipa.close()


# Main
parser = argparse.ArgumentParser(description='Split input fasta file based on\
                                              existence of blast hits in input\
                                              blast table file')
parser.add_argument('-p', action='store_true', help='switch for parameter file\
                                                     output')
parser.add_argument('-f', type=str, help='fastafile')
parser.add_argument('-b', type=str, help='blastfile')
args = parser.parse_args()
try:
    binf = open(args.b)
except IOError:
    sys.exit('input file error')
faLst = FastaList(args.f)
finf = faLst.faFile
spFaLst = spadesFa(finf)
blLst = blastTbl(binf)
WrFiles()
finf.close()
if args.f[-2:] == 'gz' and os.path.isfile(args.f[:-3]):
    os.remove(args.f[:-3])
binf.close()
