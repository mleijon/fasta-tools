#!/uar/bin/python3

class genbank_parse(object):
    """docstring for genbank_record."""
    mol_types = ['Genomic DNA', 'Genomic RNA', 'Precursor RNA', 'mRNA, cDNA',
                 'Ribosomal RNA', 'Transfer RNA', 'Other-Genetic', 'cRNA',
                 'Transcribed RNA', 'Transfer-messenger RNA', 'ncRNA']
    div_types = [ 'PRI', 'ROD', 'MAM', 'VRT', 'INV', 'PLN', 'BCT', 'VRL',
                  'PHG', 'SYN', 'UNA', 'EST', 'PAT', 'STS', 'GSS', 'HTG',
                  'HTC', 'ENV']
    fasta_width = 70
    def __init__(self, record):
        super(genbank_record, self).__init__()
        self.record = record
    def extract(str_org, str_1, str_2):
        l_index = str_org.find(str_1) + len(str_first)
        u_index = str_org.find(str_2)
        return str_org[l_index:u_index]
    self.LOCUS = extract(record, 'LOCUS', 'DEFINITION')
    locus_lst = record.LOCUS.split('  ').strip()
    self.LOCUS.name = locus_lst[0]
    self.LOCUS.seq_len = locus_lst[1]
    self.LOCUS.mol_type = locus_lst[2]
    self.LOCUS.topology = locus_lst[3]
    self.LOCUS.gb_div = locus_lst[4].split(' ').strip()[0]
    self.LOCUS.date = locus_lst[4].split(' ').strip()[1]
    self.DEFINITION = extract(record, 'DEFINITION', 'ACCESSION')
    self.ACCESSION = extract(record, 'ACCESSION','VERSION')
    self.SEQUENCE = ''
    counter = 0
    sequence = extract(record,'ORIGIN','//')
    for char in sequence:
        if not (char.isnumeric() or char.isspace()):
            self.SEQUENCE += char
            counter += 1
            if counter == 70:
                self.SEQUENCE += '\n'
                counter = 0
    def make_fasta(self):
        fasta = ''.join(['>', self.ACCESSION, self.DEFINITION, '\n',
                         self.SEQUENCE])
        if  not fasta.endswith('\n'):
            fasta.append('\n')
gb_file = open('genbank.gb')
record = genbank.read()
gb = gb_parse(record)
fastaout = gb.make_fasta()
print(fastaout)
