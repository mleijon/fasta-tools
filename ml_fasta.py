#Creates a list of the sequence ids from input fasta file
def rd_ids (fa_file):
    id_list = []
    for line in fa_file:
        if line.startswith('>'):
            id_list.append(line[1:].strip())
    fa_file.seek(0)
    return id_list

#creates a list of fasta test sequences from inut fasta file
def rd_seqs (fa_file):
    seq_list = []
    newseq = ''
    first = True
    for line in fa_file:
        if line.startswith('>'):
            if not first:
                seq_list.append(newseq)
            newseq = line
        else:
            first = False
            newseq +=line
    seq_list.append(newseq)
    nr_seq = len(seq_list)
    fa_file.seek(0)
    return seq_list
