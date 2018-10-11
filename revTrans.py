#!/usr/bin/python3
class revTrans(object):
    """docstring for revTrans."""

    genCode = {"F":["TTT","TTC"],"L":["TTA","TTG", "CTT", "CTC", "CTA", "CTG"],"A":["GCT", "GCC", "GCA", "GCG"],
               "R":["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],"N":["AAT", "AAC"],"D":["GAT", "GAC"],
               "C":["TGT", "TGC"],"Q":["CAA", "CAG"],"E":["GAA", "GAG"],"G":["GGU", "GGC", "GGA", "GGG"],
               "H":["CAT", "CAC"],"I":["ATT", "ATC", "ATA"],"K":["AAA", "AAG"],"M":["ATG"],
               "P":["CCU", "CCC", "CCA", "CCG"],"S":["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
               "T":["ACT", "ACC", "ACA", "ACG"],"W":["TGG"],"Y":["TAT", "TAC"],"V":["GTT", "GTC", "GTA", "GTG"],
               "*":["TAA", "TGA", "TAG"],"X":["TTT","TTC","TTA","TTG", "CTT", "CTC", "CTA", "CTG","GCT", "GCC",
                    "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG","AAT", "AAC","GAT", "GAC","TGT",
                    "TGC","CAA", "CAG","GAA", "GAG","GGU", "GGC", "GGA", "GGG","CAT", "CAC","ATT", "ATC",
                    "ATA","AAA", "AAG","ATG","CCU", "CCC", "CCA", "CCG","TCT", "TCC", "TCA", "TCG", "AGT",
                    "AGC","ACT", "ACC", "ACA", "ACG","TGG","TAT", "TAC","GTT", "GTC", "GTA", "GTG"]}

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
        if pep == "":
            self.genLst = genLst
        else:
            crGen(pep, genLst)

#main
fusionPeptide = 'GLFGAIAGFI'
gl = revTrans(fusionPeptide)

def parseSeq (seq):
    pass

def wrfafromls (seqLst):
    j = 0
    filFa = open('seqfile.fa','w')
    for seq in seqLst:
        j+=1
        filFa.write('> %d\n' % j)
        filFa.write(seq+'\n')
    filFa.close()

wrfafromls(gl.genLst)
