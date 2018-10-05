#!/usr/bin/python3

#Creates a list of dictionaries of the sequence parameters assuming SPADES (http://cab.spbu.ru/software/spades/) contigs
# NODE (the sequence number), length (length of sequence), cov (read coverage of the contig)
def rd_par_fas(fa_file):
    t = []
    seq_pars = dict()
    seq_par_list = []
    for line in fa_file:
        if line.startswith('>'):
            t = line.split('_')
            seq_pars ={t[0][1:]:int(t[1]),t[2]:int(t[3]),t[4]:float(t[5])}
            seq_par_list.append(seq_pars)
    return seq_par_list

def countlines(infile):
    with infile as f:
        for count, element in enumerate(f,1):
            pass
    return count

#Creates dictionary of blast *no hits* ids
def rd_ids_bla(bla_file):
    global nr_rmd
    no_hit = '***** No hits found ******'
    linelist = []
    no_hits = dict()
    j = 0
    for line in bla_file:
        j+=1
        if line.strip() == no_hit:
            linelist.append(j-6)
    for linenr in linelist:
        no_hits[linecache.getline(args.b,linenr).strip()[7:]] = linenr
    nr_rmd = len(linelist)
    linecache.clearcache()
    bla_file.seek(0)
    return no_hits

# Write filtered fasta file not including sequeces with no blast hits
def write_filtered_fa(fasta_in, blast_in):
    fil_fa = open('hits_'+args.f,'w')
    fil_re = open('nohits_'+args.f,'w')
    fasta_in_seqs = rd_seqs(fasta_in)
    fasta_in_seqids = rd_ids(fasta_in)
    fasta_exclude = rd_ids_bla(blast_in)
    j = 0
    for ids in fasta_in_seqids:
        if ids not in fasta_exclude:
            fil_fa.write(fasta_in_seqs[j])
        else:
            fil_re.write(fasta_in_seqs[j])
        j+=1
    fil_fa.close()
    fil_re.close()

#Check if the file appears to be created by Spades
def spades_file(fa_file):
    parameters = rd_par_fas(fa_file)
    for i in range(len(parameters)):
        if 'NODE' not in parameters[i]:
            return False
        elif 'length' not in parameters[i]:
            return False
        elif 'cov' not in parameters[i]:
            return False
    fa_file.seek(0)
    return True

#Writes a parameter file based on the ID-string created by Spades and som other parameters
def write_parfile (fa_file,fipa,inpname):
    parameters = rd_par_fas(fa_file)
    nr_of_contigs = len(parameters)
    sum_of_cov = 0
    sum_of_len = 0
    fipa.write('Parameters for: '+inpname+'\n\n')
    fipa.write('Number of contigs: %d\n' % nr_of_contigs)
    for i in range(len(parameters)):
        sum_of_cov += parameters[i]['cov']
        sum_of_len += parameters[i]['length']
    aver_length = round(sum_of_len/nr_of_contigs,0)
    aver_cov = round(sum_of_cov/nr_of_contigs,0)
    fipa.write('Average contig length: %d\n' % aver_length)
    fipa.write('Average read coverage: %d\n' % aver_cov)
    fipa.write('_________________________________________\n\n')

#Main
import linecache
import argparse
import sys
from ml_fasta import rd_ids, rd_seqs
parser = argparse.ArgumentParser(description='Split input fasta file based on existence of blast hits in input blast table file')
parser.add_argument('-p',action='store_true',help='switch for parameter file output')
parser.add_argument('-f',type = str,help='fastafile')
parser.add_argument('-b',type = str,help='blastfile')
args = parser.parse_args()
try:
    finf = open(args.f)
except:
    sys.exit('input file error')
try:
    binf = open(args.b)
except:
    sys.exit('input file error')
try:
    write_filtered_fa(finf,binf)
except:
    finf.close()
    binf.close()
    sys.exit('error writing file')
if args.p and spades_file(finf):
    if '.' in args.f:
        name = args.f[:args.f.find('.')]+'.par'
        fipa = open(name,'w')
    else:
        fipa = open(args.f+'.par','w')
    write_parfile(finf,fipa,args.f)
    fihi = open('hits_'+args.f)
    write_parfile(fihi,fipa,'hits_'+args.f)
    finf.seek(0)
    percent_remove = round(100*nr_rmd/len(rd_seqs(finf)),0)
    fipa.write('%d sequences removed (%g%%)' % (nr_rmd,percent_remove))
    fipa.close()
    fihi.close()
elif args.p:
    fipa = open('par_'+args.f,'w')
    fipa.write('No parameters, not a contigs fasta file created by spades.')
    fipa.close()
print(countlines(binf))
finf.close()
binf.close()
