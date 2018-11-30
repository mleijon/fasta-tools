#!/usr/bin/python3


class GbParse(object):
    """docstring for genbank_record."""

    mol_types1 = ['DNA', 'RNA', 'mRNA', 'cDNA', 'Other-Genetic', 'cRNA',
                  'ncRNA', 'ss-RNA', 'ds-RNA', 'ss-DNA', 'ds-DNA',
                  'tRNA', 'ds-cRNA', 'ms-DNA']
    mol_types2 = ['Genomic DNA', 'Genomic RNA', 'Precursor RNA',
                  'Ribosomal RNA', 'Transfer RNA', 'Transcribed RNA',
                  'Transfer-messenger RNA']
    div_types = ['PRI', 'ROD', 'MAM', 'VRT', 'INV', 'PLN', 'BCT', 'VRL',
                 'PHG', 'SYN', 'UNA', 'EST', 'PAT', 'STS', 'GSS', 'HTG',
                 'HTC', 'ENV']
    fasta_width = 70

    def __init__(self, record):
        super(GbParse, self).__init__()

        def extract(str_org, str_1, str_2):
            l_index = str_org.find(str_1) + len(str_1)
            u_index = str_org.find(str_2)
            return str_org[l_index:u_index]

        self.LOCUS = self.del_exblanks(extract(record, 'LOCUS', 'DEFINITION'))
        locus_lst = list(filter(lambda x: x != '', self.LOCUS.split(' ')))
        #print(locus_lst)
        self.name = locus_lst[0].strip()
        self.seq_len = locus_lst[1].strip() + ' ' + locus_lst[2].strip()
        if locus_lst[3].strip() in self.mol_types1:
            self.mol_type = locus_lst[3].strip()
            j = 0
        else:
            self.mol_type = locus_lst[3].strip() + locus_lst[4].strip()
            j = 1
        self.topology = locus_lst[4 + j].strip()
        self.gb_div = locus_lst[5 + j].strip()
        self.date = locus_lst[6 + j].strip()
        self.DEFINITION = extract(record, 'DEFINITION', 'ACCESSION').strip()
        self.DEFINITION = self.DEFINITION.replace('\n', '')
        self.DEFINITION = self.del_exblanks(self.DEFINITION)[:-1]
        if 'DBLINK' in extract(record, 'VERSION', 'KEYWORDS'):
            self.VERSION = extract(record, 'VERSION', 'DBLINK').strip()
        else:
            self.VERSION = extract(record, 'VERSION', 'KEYWORDS').strip()
        self.SEQUENCE = ''
        counter = 0
        sequence = extract(record, 'ORIGIN', '//\n').split('\n', 1)[1]
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
