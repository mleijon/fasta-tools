import sys


class FastaList():
    """docstring for fasta."""

    def __init__(self, fasta_name):  # Initialized by a filename string
        seq_list = []
        id_list = []
        newseq = ''
        self.is_fastq = True  # File ca be gzipped or fastq as well as fasta
        self.name = fasta_name
        self.is_gzip = self.name[-2:].casefold() == 'gz'
        if '.fq' in self.name.casefold():
            self.fq_ext = '.fq'
        elif '.fastq' in self.name.casefold():
            self.fq_ext = '.fastq'
        else:
            self.is_fastq = False
        if self.is_fastq:
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
        if self.is_gzip:
            faFiuz = open(self.name[:-3], 'w')
            with gzip.open(self.name, 'rt') as f:
                faFiuz.write(f.read())
            faFiuz.close()
            faFiuz = open(self.name[:-3])
            return faFiuz
        else:
            faFi = open(self.name)
            return faFi

    def fq2fa(self):
        outFa = open(self.name[:self.name.find(self.fq_ext)]+'.fa', 'w')
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
                outFa.close()
                sys.exit('Not a fastq-file')
            j += 1
        outFa.close()
        outFa = open(self.name[:self.name.find(self.fq_ext)]+'.fa')
        return outFa

    def rev(self, seqLst):
        return seqLst[::-1]

    def revComp(self, seqLst):
        pass
