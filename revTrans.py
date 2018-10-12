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

class seqVar(object):
    """docstring for seqVar."""
    def __init__(self, peptide):
        self.pepStr = peptide
    def cheStr(self):
        ctL = 0; ctR = 0
        allowed = ["A","C","D","E","F","G","H","I","K","L","M","N",
        "P","Q","R","S","T","V","W","X","Y","*","[","]"]
        for ch in self.pepStr:
            if ch not in allowed:
                return False
            elif ch == "[":
                ctL += 1
            elif ch == "]":
                ctR += 1
        if ctL != ctR:
            return False
        else:
            return True
    def crVar(self):
        genLst = []
        variants = []
        varLst = []
        frag = ""
        tmpLst = []
        degen = False
        for ch in self.pepStr:
            if ch == "]":
                degen = False
                varLst.append(variants)
                variants = []
            elif degen:
                variants.append(ch)
            else:
                if ch == "[":
                    if frag != "":
                        tmpLst.append(frag)
                        varLst.append(tmpLst)
                        frag = ""
                        tmpLst = []
                    degen = True
                else:
                    frag += ch
        if frag != "":
            tmpLst.append(frag)
            varLst.append(tmpLst)
            tmpLst = []


        for elem in varLst:
            tmpLst = genLst.copy()
            genLst = []
            for frag in elem:
                for segment in tmpLst:
                    genLst.append(segment + frag)
        return genLst


#main
fusionPeptide = 'GLFGA[GHW]IAGFI'
gl = seqVar(fusionPeptide)

print(gl.crVar())


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

#wrfafromls(gl.genLst)
