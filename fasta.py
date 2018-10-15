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
                print('Decompressing...',end='',flush=True)
            faFiuz.close()
            faFiuz = open(self.faFile[:-3])
            print('Done decompressing',end='',flush=True)
            return faFiuz
        else:
            faFi = open(self.faFile)
            return faFi
    def fq2fa(self):
        outFa = open(self.faFile[:self.faFile.find('.gz')]+'fa')
        j = 0
        for line in self.rdfi():
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
                break
                outFa.close()
        outFa.close()

    def rev(self,seqLst):
        return seqLst[::-1]
    def revComp(self,seqLst):
        pass
