class FastaList(object):
    """docstring for fasta."""
    def __init__(self, fastafile):
        seq_list = []
        id_list = []
        newseq = ''
        first = True
        for line in fastafile:
            if line.startswith('>'):
                id_list.append(line[1:].strip())
                if not first:
                    seq_list.append(newseq)
                newseq = line
            else:
                first = False
                newseq +=line
        seq_list.append(newseq)
        fastafile.seek(0)
        nr_seq = len(seq_list)
        self.seq_list = seq_list
        self.id_list = id_list
        self.nr_seq = len(id_list)
