#!/usr/bin/python3
"""General purpose module to handle a fasta file object. Converst fastq to
fasta and unzip gzipped files"""
import sys


class FastaList():
    """docstring for fasta."""

    def __init__(self, fasta_name):  # Initialized by a filename string
        seq_list = []
        id_list = []
        newseq = ''
        self.name = fasta_name
        self.is_gzip = self.name.endswith(('.gz', '.GZ'))
        if (self.is_gzip and self.name[:-3].endswith(('.fq', '.FQ', '.fastq',
                                                      '.FASTQ')) or self.name.
                endswith(('.fq', '.FQ', '.fastq', '.FASTQ'))):
            self.fa_file = self.fq2fa()
        else:
            self.fa_file = self.rdfi()
        for line in self.fa_file:
            if line.startswith('>'):
                if id_list != []:
                    seq_list.append(newseq + '\n')
                id_list.append(line[1:].strip())
                newseq = line
            else:
                newseq += line.strip()
        seq_list.append(newseq + '\n')
        self.fa_file.seek(0)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)

    def rdfi(self):
        """ Reads and return a fasta file unzipped if necessary"""
        import gzip
        if self.is_gzip:
            fa_fiuz = open(self.name[:-3], 'w')
            with gzip.open(self.name, 'rt') as in_fi:
                fa_fiuz.write(in_fi.read())
            fa_fiuz.close()
            fa_fiuz = open(self.name[:-3])
            return fa_fiuz
        fa_fi = open(self.name)
        return fa_fi

    def fq2fa(self):
        """ Reads and return a fastq file converted to fasta format"""
        outfa_name = self.name.rpartition('/')[2][:self.name.rpartition('/')
                                                  [2].find('.')] + '.fa'
        outfa_name = self.name.rpartition('/')[0] + '/' + outfa_name
        out_fa = open(outfa_name, 'w')
        j = 0
        for line in self.rdfi():
            if line[0] == '@' and j % 4 == 0:
                out_fa.write(line.replace('@', '>', 1))
            elif j % 4 == 1:
                out_fa.write(line)
            elif line[0] == '+' and j % 4 == 2:
                pass
            elif j % 4 == 3:
                pass
            else:
                out_fa.close()
                sys.exit('Not a fastq-file')
            j += 1
        out_fa.close()
        out_fa = open(outfa_name)
        return out_fa

    def rev(self):
        """Yields the reverse of a nucleotide sequence"""
        return self.seq_list[::-1]

    def rev_comp(self, seq_lst):
        """Yields the reverse complement of a nucleotide sequence"""
        pass
