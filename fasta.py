#!/usr/bin/python3
"""General purpose module to handle a fasta file object. Converts fastq to
fasta and unzip gzipped files"""
import sys


class FastaList:
    """From a fasta text file (fastq or fasta) Creates a list with the id:s of
    the fasta sequences (id_list) and a list of the sequences with the newlines
    removed (seq_list), suitable for searching. The sequenes list elements are
    strings with a newline separating the id and the nucleotides sequence."""

    def __init__(self, fasta_name):  # Initialized by a filename string
        self.seq_list = []
        self.id_list = []
        newseq = ''
        self.name = fasta_name
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
        # Handles input with full path
        basename = self.name.rpartition('/')[2]
        outfa_name = basename[:basename.find('.')] + '.fa'
        outfa_name = self.name.rpartition('/')[0] + '/' + outfa_name
        out_fa = open(outfa_name, 'w')
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
        out_fa = open(outfa_name)
        return out_fa

    def rev(self):
        """Returns the reverse of a nucleotide sequence"""
        return self.seq_list[::-1]

    def rev_comp(self):
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
            fasta_id = fasta.split('\n')[0]
            fasta_seq = fasta.split('\n')[1]
            fasta_seq_rev = fasta_seq[::-1]
            new_fasta = fasta_id + '_RC\n'
            for nucl in fasta_seq_rev:
                new_fasta += comp(nucl)
            new_fasta += '\n'
            seq_list_rc.append(new_fasta)
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
