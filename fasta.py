#!/usr/bin/python3
"""General purpose module to handle a fasta file object. Converts fastq to
fasta and unzip gzipped files"""
import sys
import os


class FastaList:
    """From a fasta text file (fastq or fasta) Creates a list with the id:s of
    the fasta sequences (id_list) and a list of the sequences with the newlines
    removed (seq_list), suitable for searching. The sequenes list elements are
    strings with a newline separating the id and the nucleotides sequence.
    id_list elements will end without \n while seq_list elements ends with \n.
    """

    def __init__(self, fasta_name):  # Initialized by a filename string
        self.seq_list = []
        self.id_list = []
        newseq = ''
        self.name = str(os.path.abspath(fasta_name))
        self.is_gzip = self.name.endswith(('.gz', '.GZ'))
        if (self.is_gzip and self.name[:-3].endswith(('.fq', '.FQ', '.fastq',
                                                      '.FASTQ')) or self.name.
                endswith(('.fq', '.FQ', '.fastq', '.FASTQ'))):
            self.fa_file = self.fq2fa()  # Handle file as fastq
        else:
            self.fa_file = self.rdfi()  # Handle file as fasta
        for line in self.fa_file:
            if line.startswith('>'):
                if self.id_list:
                    self.seq_list.append(newseq + '\n')
                self.id_list.append(line[1:].strip())
                newseq = line
            else:
                newseq += line.strip()
        self.seq_list.append(newseq + '\n')
        self.fa_file.seek(0)
        self.nr_seq = len(self.id_list)

    def rdfi(self):
        """ Reads and return a fasta file object unzipped if necessary"""
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
        out_fa = open(self.name[:self.name.find('.')] + '.fa', 'w')
        for line in enumerate(self.rdfi()):
            if line[1][0] == '@' and line[0] % 4 == 0:
                out_fa.write(line[1].replace('@', '>', 1))  # seqid line
            elif line[0] % 4 == 1:
                out_fa.write(line[1])  # sequence line
            elif line[1][0] == '+' and line[0] % 4 == 2:
                pass
            elif line[0] % 4 == 3:
                pass
            else:
                out_fa.close()
                sys.exit('Not a fastq-file')
        out_fa.close()
        out_fa = open(out_fa.name)
        return out_fa

    def seq_list_rev(self):
        """Returns the reverse of a nucleotide sequence"""
        rev_seq_list = []
        for item in self.seq_list:
            seqid = item.split('\n')[0] + '_rev\n'
            seq = item.split('\n')[1][::-1] + '\n'
            rev_seq_list.append(seqid + seq)
        return rev_seq_list

    def seq_list_revc(self):
        """Creates a list the reverse complement of a nucleotide sequences in
        the list seq_list"""
        seq_list_rc = []

        def comp(nt):
            return {
                'A': 'T',
                'T': 'A',
                'G': 'C',
                'C': 'G',
                'Y': 'R',
                'R': 'Y',
                'K': 'M',
                'M': 'K',
                'S': 'S',
                'W': 'W',
                'N': 'N',
                'B': 'V',
                'V': 'B',
                'D': 'H',
                'H': 'D'
            }[nt.upper()]
        for fasta in self.seq_list:
            seqid = fasta.split('\n')[0] + '_rc\n'
            seq = fasta.split('\n')[1]
            seq_rc = ''
            seq_rev = seq[::-1]
            for nucl in seq_rev:
                seq_rc += comp(nucl)
            seq_rc += '\n'
            seq_list_rc.append(seqid + seq_rc)
        return seq_list_rc

    def divide(self, divisor):
        tmp_list = self.seq_list.copy()
        seq_list_div = []
        item_per_lst = self.nr_seq//divisor + (self.nr_seq % divisor) //\
            divisor
        for i in range(divisor - 1):
            seq_list_div.append(tmp_list[0:item_per_lst])
            del tmp_list[0:item_per_lst]
        seq_list_div.append(tmp_list)
        return seq_list_div


if __name__ == "__main__":
    testfi = FastaList('ml.fa')
    for j in range(1):
        print(testfi.seq_list[j])
    for j in range(1):
        print(testfi.seq_list_rev()[j])
    for j in range(1):
        print(testfi.seq_list_revc()[j])
