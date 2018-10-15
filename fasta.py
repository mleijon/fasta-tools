class FastaList(object):
    """docstring for fasta."""

    def __init__(self, faFile):
        seq_list = []
        id_list = []
        newseq = ''
        self.faFile = faFile
        first = True
        for line in rdfi():
            if line.startswith('>'):
                id_list.append(line[1:].strip())
                if not first:
                    seq_list.append(newseq)
                newseq = line
            else:
                first = False
                newseq +=line
        seq_list.append(newseq)
        rdfi().seek(0)
        nr_seq = len(seq_list)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)
    def rdfi(self):
        import gzip; import tempfile

        if self.faFile[-2:] == 'gz':
            ftemp = tempfile.TemporaryFile()
            with gzip.open(self.faFile,'rt') as f:
                ftemp.write(f.read())
            return ftemp
        else:
            fastafile = open(self.faFile)
            return fastafile
    def rev(self,seqLst):
        return seqLst[::-1]
    def revComp(self,seqLst):
        pass
