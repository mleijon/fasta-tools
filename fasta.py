class FastaList(object):
    """docstring for fasta."""

    def __init__(self, faFile):
        seq_list = []
        id_list = []
        newseq = ''
        self.faName = faFile
        if self.faName[-2:] == 'gz':
            self.gz = True
        if self.gz:
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
                newseq +=line
        seq_list.append(newseq)
        self.rdfi().seek(0)
        nr_seq = len(seq_list)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)

    def rdfi(self):
        import gzip
        if self.gz:
            faFiuz = open(self.faName[:-3],'w')
            with gzip.open(self.faName,'rt') as f:
                faFiuz.write(f.read())
                print('Decompressing...',end='',flush=True)
            faFiuz.close()
            faFiuz = open(self.faName[:-3])
            print('Done decompressing',end='',flush=True)
            return faFiuz
        else:
            faFi = open(self.faName)
            return faFi
    def fq2fa(self):
        outFa = open(self.faName[:self.faName.find('.gz')]+'fa')
        j = 0
        for line in self.rdfi():
            print('converting...',end='',flush=True)
            if line[0] == '@' and j % 4 == 0:
                outFa.write(line.replace('@','>',1))
            elif j % 4 == 1:
                outFa.write(line)
            elif line[0] == '+' and j % 4 == 2:
                pass
            elif j % 4 == 3:
                pass
            else
                print('Not a fastq-file')
                return
                outFa.close()
        print('Done converting...',end='',flush=True)
        outFa.close()
        return outFa

    def rev(self,seqLst):
        return seqLst[::-1]
    def revComp(self,seqLst):
        pass
