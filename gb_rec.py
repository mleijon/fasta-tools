#!/usr/bin/python3


class genbank_parse(object):
    """docstring for genbank_record."""
    mol_types = ['Genomic DNA', 'Genomic RNA', 'Precursor RNA', 'mRNA, cDNA',
                 'Ribosomal RNA', 'Transfer RNA', 'Other-Genetic', 'cRNA',
                 'Transcribed RNA', 'Transfer-messenger RNA', 'ncRNA']
    div_types = ['PRI', 'ROD', 'MAM', 'VRT', 'INV', 'PLN', 'BCT', 'VRL',
                 'PHG', 'SYN', 'UNA', 'EST', 'PAT', 'STS', 'GSS', 'HTG',
                 'HTC', 'ENV']
    fasta_width = 71

    def __init__(self, record):
        super(genbank_parse, self).__init__()

        def extract(str_org, str_1, str_2):
            l_index = str_org.find(str_1) + len(str_1)
            u_index = str_org.find(str_2)
            return str_org[l_index:u_index]

        def del_exblanks():
            old_len = len(self.DEFINITION)
            self.DEFINITION = self.DEFINITION.replace('  ', ' ')
            new_len = len(self.DEFINITION)
            if not(new_len == old_len):
                del_exblanks()
                return self.DEFINITION

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
        self.DEFINITION = del_exblanks()
        self.ACCESSION = extract(record, 'ACCESSION', 'VERSION').strip()
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

    def make_fasta(self):
        fasta = ''.join(['>', self.ACCESSION, ' ', self.DEFINITION, '\n',
                         self.SEQUENCE])
        if not fasta.endswith('\n'):
            fasta += '\n'
        return fasta


gb_file = open('sequence.gb')
record = gb_file.read()
gb = genbank_parse(record)
fastaout = gb.make_fasta()
print(fastaout)
