#!/usr/bin/python3

#Creates a list of dictionaries of the sequence parameters assuming SPADES (http://cab.spbu.ru/software/spades/) contigs
# NODE (the sequence number), length (length of sequence), cov (read coverage of the contig)
from fasta import FastaList
class spadesFa(FastaList):
    """docstring for ."""
    def __init__(self, spades_file):
        super(spadesFa,self).__init__(spades_file)
        sumCov = 0
        sumLen = 0
        t = []
        seqPar = dict()
        self.seqParLst = []
        for line in spades_file:
            if line.startswith('>'):
                t = line.split('_')
                seqPar ={t[0][1:]:int(t[1]),t[2]:int(t[3]),t[4]:float(t[5])}
                self.seqParLst.append(seqPar)
        self.nr_of_contigs = len(self.seqParLst)
        for i in range(len(self.seqParLst)):
            sumCov += self.seqParLst[i]['cov']
            sumLen += self.seqParLst[i]['length']
        self.aver_length = round(sumLen/self.nr_of_contigs,0)
        self.aver_cov = round(sumCov/self.nr_of_contigs,0)
        self.is_spadesFa = True
        for i in range(len(self.seqParLst)):
            if 'NODE' not in self.seqParLst[i]:
                self.is_spadesFa = False
            elif 'length' not in self.seqParLst[i]:
                self.is_spadesFa = False
            elif 'cov' not in self.seqParLst[i]:
                self.is_spadesFa = False

#Writes a parameter file based on the ID-string created by Spades and some other parameters
    def wrPar2File(self,parfile,inpname):
            parfile.write('Parameters for: '+inpname+'\n\n')
            parfile.write('Number of contigs: %d\n' % self.nr_of_contigs)
            parfile.write('Average contig length: %d\n' % self.aver_length)
            parfile.write('Average contig coverage: %d\n' % self.aver_cov)
            parfile.write('_________________________________________\n\n')

#Creates dictionary of blast *no hits* ids
def rd_ids_bla(bla_file):
    global nr_rmd       #nr of contigs removed
    no_hit = '***** No hits found ******'
    linelist = []
    no_hits = dict()
    j = 0
    for line in bla_file:
        j+=1
        if line.strip() == no_hit:
            no_hits[linecache.getline(args.b,j-6).strip()[7:]] = j-6
    nr_rmd = len(no_hits)
    linecache.clearcache()
    bla_file.seek(0)
    return no_hits

# Write filtered fasta file not including sequeces with no blast hits
def wrFiltFa(blast_in):
    fil_fa = open('hits_'+args.f,'w')
    fil_re = open('nohits_'+args.f,'w')
    fasta_exclude = rd_ids_bla(blast_in)
    j = 0
    for ids in faLst.id_list:
        if ids not in fasta_exclude:
            fil_fa.write(faLst.seq_list[j])
        else:
            fil_re.write(faLst.seq_list[j])
        j+=1
    fil_fa.close()
    fil_re.close()

#Main
import linecache
import argparse
import sys
from fasta import FastaList
parser = argparse.ArgumentParser(description='Split input fasta file based on existence of blast hits in input blast table file')
parser.add_argument('-p',action='store_true',help='switch for parameter file output')
parser.add_argument('-f',type = str,help='fastafile')
parser.add_argument('-b',type = str,help='blastfile')
args = parser.parse_args()
try:
    finf = open(args.f)
    binf = open(args.b)
except:
    sys.exit('input file error')
faLst = spadesFa(finf)
try:
    wrFiltFa(binf)
except:
    binf.close()
    sys.exit('error writing file')
if args.p and faLst.is_spadesFa:
    if '.' in args.f:
        name = args.f[:args.f.find('.')]+'.par'
        fipa = open(name,'w')
    else:
        fipa = open(args.f+'.par','w')
    faLst.wrPar2File(fipa,args.f)
    fihi = open('hits_'+args.f)
    faHitLst = spadesFa(fihi)
    faHitLst.wrPar2File(fipa,'hits_'+args.f)
    percent_remove = round(100*nr_rmd/len(faLst.seq_list),0)
    fipa.write('%d sequences removed (%g%%)' % (nr_rmd,percent_remove))
    fipa.close()
    fihi.close()
elif args.p:
    fipa = open('par_'+args.f,'w')
    fipa.write('No parameters, not a contigs fasta file created by spades.')
    fipa.close()
finf.close()
binf.close()
