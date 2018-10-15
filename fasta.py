class FastaList(object):
    """docstring for fasta."""

    def __init__(self, faFile):
        seq_list = []
        id_list = []
        newseq = ''
        self.faFile = faFile
        first = True
        for line in self.rdfi():
            if line.startswith('>'):
                id_list.append(line[1:].strip())
                if not first:
                    seq_list.append(newseq)
                newseq = line
            else:
                first = False
                newseq +=line
        seq_list.append(newseq)
        self.rdfi().seek(0)
        nr_seq = len(seq_list)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)
    def rdfi(self):
        import gzip

        if self.faFile[-2:] == 'gz':
            faFiuz = open(self.faFile[:-3],'w')
            with gzip.open(self.faFile,'rt') as f:
                faFiuz.write(f.read())
            faFiuz.close()
            faFiuz = open(self.faFile[:-3])
            return faFiuz
        else:
            faFi = open(self.faFile)
            return faFi

    def rev(self,seqLst):
        return seqLst[::-1]
    def revComp(self,seqLst):
        pass
