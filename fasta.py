class FastaList(object):
    """docstring for fasta."""
    def __init__(self, fa_file):
        seq_list = []
        id_list = []
        newseq = ''
        first = True
        for line in fa_file:
            if line.startswith('>'):
                id_list.append(line[1:].strip())
                if not first:
                    seq_list.append(newseq)
                newseq = line
            else:
                first = False
                newseq +=line
        seq_list.append(newseq)
        nr_seq = len(seq_list)
        fa_file.seek(0)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)
