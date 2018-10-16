class FastaList(object):
    """docstring for fasta."""

    def __init__(self, faFile):
        seq_list = []; id_list = []; newseq = ''
        self.gz = False; self.fq = True
        self.faName = faFile
        if self.faName[-2:] in ['gz','GZ']:
            self.gz = True
        if '.fq' in self.faName:
            self.fqExt = '.fq'
        elif '.FQ' in self.faName:
            self.fqExt = '.FQ'
        elif '.fastq' in self.faName:
            self.fqExt = '.fastq'
        elif '.FASTQ' in self.faName:
            self.fqExt = '.FASTQ'
        else:
            self.fq = False
        if '.gz' in self.faName:
            self.gzExt = 'gz'
        elif '.GZ' in self.faName:
            self.gzExt = 'GZ'
        else:
            self.gz = False

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
                newseq +=line
        seq_list.append(newseq)
        self.faFile.seek(0)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)

    def rdfi(self):
        import gzip
        if self.gz:
            faFiuz = open(self.faName[:-3],'w')
            with gzip.open(self.faName,'rt') as f:
                faFiuz.write(f.read())
            faFiuz.close()
            faFiuz = open(self.faName[:-3])
            return faFiuz
        else:
            faFi = open(self.faName)
            return faFi
    def fq2fa(self):
        if self.gz:
            outFa = open(self.faName[:self.faName.find(self.fqExt + '.gz')]+'.fa','w')
        else:
            outFa = open(self.faName[:self.faName.find(self.fqExt)]+'.fa','w')
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
            else:
                print('Not a fastq-file')
                return
                outFa.close()
            j+=1
        outFa.close()
        if self.gz:
            outFa = open(self.faName[:self.faName.find(self.fqExt + '.gz')]+'.fa')
        else:
            outFa = open(self.faName[:self.faName.find(self.fqExt)]+'.fa')
        return outFa

    def rev(self,seqLst):
        return seqLst[::-1]
    def revComp(self,seqLst):
        pass
