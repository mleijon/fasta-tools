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
        super(spadesFa, self).__init__(ARGS.f)
        sum_cov = 0
        sum_len = 0
        t = []
        seq_par = dict()
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
                seq_par = {t[0][1:]: int(t[1]), t[2]: int(t[3]), t[4]:
                           float(t[5])}  # parses the id strings of fa-file
                self.seqParLst.append(seq_par)
        if self.is_spadesFa:
            self.nr_contig = len(self.seqParLst)
            for i in range(len(self.seqParLst)):
                sum_cov += self.seqParLst[i]['cov']
                sum_len += self.seqParLst[i]['length']
            self.av_len = round(sum_len/self.nr_contig, 0)
            self.av_cov = round(sum_cov/self.nr_contig, 0)

# Writes a parameter file based on the ID-string created by Spades and some
# other parameters
    def wrPar2File(self, parfile, inpname):
        parfile.write('Parameters for: '+inpname+'\n\n')
        parfile.write('Number of contigs: %d\n' % self.nr_contig)
        parfile.write('Average contig length: %d\n' % self.av_len)
        parfile.write('Average contig coverage: %d\n' % self.av_cov)
        parfile.write('_________________________________________\n\n')


class blastTbl():
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
                no_hits[linecache.getline(ARGS.b, j-6).strip()[7:]] = j-6
        linecache.clearcache()
        return no_hits

# Write filtered fasta file not including sequeces with no blast hits and a
# fasta file including exlusively these contigs. A parameter file is also
# written using the FastaList class


def WrFiles():
    if ARGS.f[-2:] in ['gz', 'GZ']:
        name_part = ARGS.f[:-3].rpartition('/')
        fahits = name_part[0] + '/' + 'hits_' + name_part[2]
        fanohits = fahits.replace('hit', 'nohit')
    else:
        name_part = ARGS.f.rpartition('/')
        fahits = name_part[0] + '/' + 'hits_' + name_part[2]
        fanohits = fahits.replace('hit', 'nohit')
    if fahits[-2:] in ['fq', 'FQ', 'fa', 'FA']:
        fahits = fahits[:-3]
        fanohits = fanohits[:-3]
    elif fahits[-5:] in ['fastq', 'FASTQ', 'fasta', 'FASTA']:
        fahits = fahits[:-6]
        fanohits = fanohits[:-6]
    fahits = fahits + '.fa'
    fanohits = fanohits + '.fa'
    try:
        FiHit = open(fahits, 'w')
        FiNoHit = open(fanohits, 'w')
    except IOError:
        sys.exit('error writing file')
    j = 0
    for ids in SPAFA_LST.id_list:
        if ids not in blLst.noHits:
            FiHit.write(SPAFA_LST.seq_list[j])
        else:
            FiNoHit.write(SPAFA_LST.seq_list[j])
        j += 1
    FiHit.close()
    FiNoHit.close()
    if ARGS.p and SPAFA_LST.is_spadesFa:
        fipa = open(fahits[:-2] + 'par', 'w')
        SPAFA_LST.wrPar2File(fipa, ARGS.f)
        fihi = open(fahits)
        faHitLst = spadesFa(fihi)
        faHitLst.wrPar2File(fipa, fahits)
        percent_remove = round(100*blLst.nrNoHits/len(SPAFA_LST.seq_list), 0)
        fipa.write('%d sequences removed (%g%%)\n' % (blLst.nrNoHits,
                                                      percent_remove))
        fipa.close()
        fihi.close()
    elif ARGS.p:
        name_par = fahits[:-3].rpartition('/')
        fipa = open(name_par[0] + '/' + 'nopar_' + name_par[2], 'w')
        fipa.write('No parameters, not contigs fasta file created by spades.')
        fipa.close()


# Main
PARSER = argparse.ArgumentParser(description='Split input fasta file based on\
                                              existence of blast hits in input\
                                              blast table file')
PARSER.add_argument('-p', action='store_true', help='switch for parameter file\
                                                     output')
PARSER.add_argument('-f', type=str, help='fastafile', required=True)
PARSER.add_argument('-b', type=str, help='blastfile', required=True)
ARGS = PARSER.parse_args()
try:
    BL_IN = open(ARGS.b)
except IOError:
    sys.exit('input file error')
faLst = FastaList(ARGS.f)
FA_IN = faLst.fa_file
SPAFA_LST = spadesFa(FA_IN)
blLst = blastTbl(BL_IN)
WrFiles()
FA_IN.close()
if ARGS.f[-2:] == 'gz' and os.path.isfile(ARGS.f[:-3]):
    os.remove(ARGS.f[:-3])
BL_IN.close()
