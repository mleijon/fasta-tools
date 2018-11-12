#!/usr/bin/python3
"""Derive all possible codons for a certain protein sequence."""
import argparse


class RevTrans():
    """docstring for RevTrans."""

    genCode = {"F": ["TTT", "TTC"],
               "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
               "A": ["GCT", "GCC", "GCA", "GCG"],
               "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
               "N": ["AAT", "AAC"],
               "D": ["GAT", "GAC"],
               "C": ["TGT", "TGC"],
               "Q": ["CAA", "CAG"],
               "E": ["GAA", "GAG"],
               "G": ["GGT", "GGC", "GGA", "GGG"],
               "H": ["CAT", "CAC"],
               "I": ["ATT", "ATC", "ATA"],
               "K": ["AAA", "AAG"],
               "M": ["ATG"],
               "P": ["CCT", "CCC", "CCA", "CCG"],
               "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
               "T": ["ACT", "ACC", "ACA", "ACG"],
               "W": ["TGG"],
               "Y": ["TAT", "TAC"],
               "V": ["GTT", "GTC", "GTA", "GTG"],
               "*": ["TAA", "TGA", "TAG"],
               "X": ["TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
                     "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG",
                     "AGA", "AGG", "AAT", "AAC", "GAT", "GAC", "TGT", "TGC",
                     "CAA", "CAG", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG",
                     "CAT", "CAC", "ATT", "ATC", "ATA", "AAA", "AAG", "ATG",
                     "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "TCA", "TCG",
                     "AGT", "AGC", "ACT", "ACC", "ACA", "ACG", "TGG", "TAT",
                     "TAC", "GTT", "GTC", "GTA", "GTG"]}

    def __init__(self, peptide):
        self.peptide = peptide
        gen_lst = []
        all_gen_lst = []
        self.valid = self.valid_str()
        if self.valid:
            for sequence in self.create_var():
                gen_lst = []
                for codon in self.genCode[sequence[0]]:
                    gen_lst.append(codon)
                pep = sequence[1:]

                def cr_genes(pep, gen_lst, all_gen_lst):
                    tmp_lst = gen_lst.copy()
                    gen_lst = []
                    for gene in tmp_lst:
                        for codon in self.genCode[pep[0]]:
                            gen_lst.append(gene + codon)
                    pep = pep[1:]
                    if pep == "":
                        all_gen_lst += gen_lst
                        return
                    cr_genes(pep, gen_lst, all_gen_lst)
                if pep != "":
                    cr_genes(pep, gen_lst, all_gen_lst)
            self.gen_lst = all_gen_lst

    def valid_str(self):
        """Control the validity of the peptide string."""
        count_left = 0
        count_right = 0
        left = False
        allowed = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "X",
                   "Y", "*", "[", "]"]
        for char in self.peptide:
            if char not in allowed:
                return False
            if char == "[":
                if left:
                    return False
                left = True
                count_left += 1
            if char == "]":
                if not left:
                    return False
                left = False
                count_right += 1
        if count_left != count_right:
            return False
        return True

    def create_var(self):
        """Creates gene variants for a given peptide seq."""
        variant = []
        var_lst = []
        tmp_lst = []
        frag = ""
        multiple = False
        for char in self.peptide:
            if char == "]":
                multiple = False
                if variant != []:
                    var_lst.append(variant)
                variant = []
            elif multiple:
                variant.append(char)
            elif char == "[":
                if frag != "":
                    var_lst.append([frag])
                    frag = ""
                multiple = True
            else:
                frag += char
        if frag != "":
            var_lst.append([frag])
        gen_lst = var_lst.pop(0)
        for elem in var_lst:
            tmp_lst = gen_lst.copy()
            gen_lst = []
            for frag in elem:
                for segment in tmp_lst:
                    gen_lst.append(segment + frag)
        return gen_lst


# main

PARSER = argparse.ArgumentParser(description='Create the genes for a peptide\
                                 taking into account degenercy of the genetic\
                                 code')
PARSER.add_argument('-p', type=str, help='amino acid sequence in single letter\
                    format', default='G[IL][FL]GA[ILK][AS]GF[IL]')
# Default = AIV fusion pep
ARGS = PARSER.parse_args()
GENE_LST = RevTrans(ARGS.p)


def wr_fasta(seq_lst):
    """ Write created gene sequences to fasta file"""
    j = 0
    fa_fi = open('seqfile.fa', 'w')
    for seq in seq_lst:
        j += 1
        fa_fi.write('> %d\n' % j)
        fa_fi.write(seq + '\n')
    fa_fi.close()


if GENE_LST.valid:
    wr_fasta(GENE_LST.gen_lst)
else:
    print("Invalid input sequence.")
