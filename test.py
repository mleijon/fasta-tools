#!/usr/bin/python3
from fasta import fasta
finf = open("test.fasta")
testfasta = fasta(finf)
print(testfasta.nr_seq)
