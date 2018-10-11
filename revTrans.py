#!/usr/bin/python3
class revTrans(object):
    """docstring for revTrans."""

    genCode = {"F":["TTT","TTC"],"L":["TTA","TTG", "CTT", "CTC", "CTA", "CTG"],"A":["GCT", "GCC", "GCA", "GCG"],
               "R":["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],"N":["AAT", "AAC"],"D":["GAT", "GAC"],
               "C":["TGT", "TGC"],"Q":["CAA", "CAG"],"E":["GAA", "GAG"],"G":["GGU", "GGC", "GGA", "GGG"],
               "H":["CAT", "CAC"],"I":["ATT", "ATC", "ATA"],"K":["AAA", "AAG"],"M":["ATG"],
               "P":["CCU", "CCC", "CCA", "CCG"],"S":["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
               "T":["ACT", "ACC", "ACA", "ACG"],"W":["TGG"],"Y":["TAT", "TAC"],"V":["GTT", "GTC", "GTA", "GTG"],
               "STOP":["TAA", "TGA", "TAG"]}

    def __init__(self, peptide):
        super(revTrans, self).__init__()
        genLst =[]
        self.peptide = peptide
        for codon in self.genCode[self.peptide[0]]:
            genLst.append(codon)
        pep = self.peptide[1:]

        def crGen(pep, genLst):
            tempLst = genLst.copy()
            genLst = []
            for gene in tempLst:
                for codon in self.genCode[pep[0]]:
                    genLst.append(gene + codon)
            pep = pep[1:]
            if pep == "":
                self.genLst = genLst
                return
            else:
                crGen(pep, genLst)
        crGen(pep, genLst)

#main
test = revTrans('MLRRGF')
print(test.genLst)
