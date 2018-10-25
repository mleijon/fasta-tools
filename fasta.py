class FastaList(object):
    """docstring for fasta."""

    def __init__(self, faFileName):  # Initialized by a filename string
        seq_list = []
        id_list = []
        newseq = ''
        self.fq = True  # File ca be gzipped or fastq as well as fasta
        self.Name = faFileName
        self.gz = self.Name[-2:].casefold() == 'gz':
        if '.fq' in self.Name.casefold():
            self.fqExt = '.fq'
        elif '.fastq' in self.Name.casefold():
            self.fqExt = '.fastq'
        else:
            self.fq = False
        if self.fq:
            self.faFile = self.fq2fa()
        else:
            self.faFile = self.rdfi()
        first = True
        for line in self.faFile:
            if line.startswith('>'):
                id_list.append(line[1:].strip())
                if not first:
                    seq_list.append(newseq)
                newseq = line
            else:
                first = False
                newseq += line
        seq_list.append(newseq)
        self.faFile.seek(0)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)

    def rdfi(self):
        import gzip
        if self.gz:
            faFiuz = open(self.Name[:-3], 'w')
            with gzip.open(self.Name, 'rt') as f:
                faFiuz.write(f.read())
            faFiuz.close()
            faFiuz = open(self.Name[:-3])
            return faFiuz
        else:
            faFi = open(self.Name)
            return faFi

    def fq2fa(self):
        outFa = open(self.Name[:self.Name.find(self.fqExt)]+'.fa', 'w')
        j = 0
        for line in self.rdfi():
            if line[0] == '@' and j % 4 == 0:
                outFa.write(line.replace('@', '>', 1))
            elif j % 4 == 1:
                outFa.write(line)
            elif line[0] == '+' and j % 4 == 2:
                pass
            elif j % 4 == 3:
                pass
            else:
                print('Not a fastq-file')
                return
                outFa.close()
            j += 1
        outFa.close()
        outFa = open(self.Name[:self.Name.find(self.fqExt)]+'.fa')
        return outFa

    def rev(self, seqLst):
        return seqLst[::-1]

    def revComp(self, seqLst):
        pass
