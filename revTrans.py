#!/usr/bin/python3
class revTrans(object):
    """docstring for revTrans."""

    genCode = {"F":["TTT","TTC"],"L":["TTA","TTG", "CTT", "CTC", "CTA", "CTG"],["GCT", "GCC", "GCA", "GCG"],
               "R":["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],"N":["AAT", "AAC"],"D":["GAT", "GAC"],
               "C":["TGT", "TGC"],"Q":["CAA", "CAG"],"E":["GAA", "GAG"],"G":["GGU", "GGC", "GGA", "GGG"],
               "H":["CAT", "CAC"],"I":["ATT", "ATC", "ATA"],"K":["AAA", "AAG"],"M":["ATG"],
               "P":["CCU", "CCC", "CCA", "CCG"],"S":["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
               "T":["ACT", "ACC", "ACA", "ACG"],"W":["TGG"],"Y":["TAT", "TAC"],"V":["GTT", "GTC", "GTA", "GTG"],
               "STOP":["TAA", "TGA", "TAG"]}

    def __init__(self, peptide):
        super(revTrans, self).__init__()
    genLst = []
    def crGen(aaPos,peptide,seqLst):
        for aa in peptide:
            for codons in geneCode[aa]:
                
