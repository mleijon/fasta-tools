#!/usr/bin/python3


class GbParse(object):
    """docstring for genbank_record."""
    # The following currently not used
    # mol_types = ['Genomic DNA', 'Genomic RNA', 'Precursor RNA', 'mRNA, cDNA',
    #              'Ribosomal RNA', 'Transfer RNA', 'Other-Genetic', 'cRNA',
    #              'Transcribed RNA', 'Transfer-messenger RNA', 'ncRNA']
    # div_types = ['PRI', 'ROD', 'MAM', 'VRT', 'INV', 'PLN', 'BCT', 'VRL',
    #              'PHG', 'SYN', 'UNA', 'EST', 'PAT', 'STS', 'GSS', 'HTG',
    #              'HTC', 'ENV']
    fasta_width = 70

    def __init__(self, record):
        super(GbParse, self).__init__()

        def extract(str_org, str_1, str_2):
            l_index = str_org.find(str_1) + len(str_1)
            u_index = str_org.find(str_2)
            return str_org[l_index:u_index]

        self.LOCUS = extract(record, 'LOCUS', 'DEFINITION')
        locus_lst = list(filter(lambda x: x != '', self.LOCUS.split('  ')))
        self.name = locus_lst[0].strip()
        self.seq_len = locus_lst[1].strip()
        self.mol_type = locus_lst[2].strip()
        self.topology = locus_lst[3].strip()
        self.gb_div = locus_lst[4].split(' ')[0].strip()
        self.date = locus_lst[4].split(' ')[1].strip()
        self.DEFINITION = extract(record, 'DEFINITION', 'ACCESSION').strip()
        self.DEFINITION = self.DEFINITION.replace('\n', '')
        self.DEFINITION = self.del_exblanks(self.DEFINITION)[:-1]
        self.VERSION = extract(record, 'VERSION', 'KEYWORD').strip()
        self.SEQUENCE = ''
        counter = 0
        sequence = extract(record, 'ORIGIN', '//')
        for char in sequence:
            if not (char.isnumeric() or char.isspace()):
                self.SEQUENCE += char.upper()
                counter += 1
                if counter == self.fasta_width:
                    self.SEQUENCE += '\n'
                    counter = 0

    def del_exblanks(self, x_str):
        old_len = len(x_str)
        x_str = x_str.replace('  ', ' ')
        new_len = len(x_str)
        while old_len > new_len:
            old_len = new_len
            x_str = x_str.replace('  ', ' ')
            new_len = len(x_str)
        return x_str

    def make_fasta(self):
        fasta = ''.join(['>', self.VERSION, ' ', self.DEFINITION, '\n',
                         self.SEQUENCE])
        if not fasta.endswith('\n'):
            fasta += '\n'
        return fasta
