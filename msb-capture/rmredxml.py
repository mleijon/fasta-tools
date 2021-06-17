#!/usr/bin/python
if __name__ == "__main__":
    import os
    import argparse
    import xml.etree.ElementTree as Et

    MAX_DEG = 10
    DEG_NUC = ['w', 's', 'm', 'k', 'r', 'y', 'b', 'd', 'h', 'v', 'n']
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-i', type=str, help='input xml-file', required=True)
    PARSER.add_argument('-o', type=str, help='output xml-file', required=True)
    ARGS = PARSER.parse_args()
    TREE = Et.parse(ARGS.i)
    ROOT = TREE.getroot()
    source_seq = {}
    source_seq_tmp = {'g'}
    nr_degseq_rm = 0
    nr_degnuc_rm = 0
    nr_no_seq = 0
    x = 0
    y = 0
    for sequence in ROOT.iter('INSDSeq'):
        y += 1
    print(y)
    for sequence in ROOT.iter('INSDSeq'):
        if sequence.find('INSDSeq_sequence') is None:
            ROOT.remove(sequence)
            nr_no_seq += 1
    for sequence in ROOT.iter('INSDSeq'):
        x += 1
    print(x)
    nr_no_seq = 0
    for sequence in ROOT.iter('INSDSeq'):
        if sequence.find('INSDSeq_sequence') is None:
            ROOT.remove(sequence)
            nr_no_seq += 1
    print(nr_no_seq)
    exit()
    for sequence in ROOT.iter('INSDSeq'):
        x += 1
        print(x)
        if sum(map(sequence.find('INSDSeq_sequence').text.count,
                   DEG_NUC)) >= MAX_DEG:
            ROOT.remove(sequence)
            nr_degnuc_rm += 1
    for sequence in ROOT.iter('INSDSeq'):
        source_seq = source_seq_tmp.copy()
        for seq in source_seq:
            if sequence.find('INSDSeq_sequence').text in seq:
                ROOT.remove(sequence)
                nr_degseq_rm += 1
                break
            elif seq in sequence.find('INSDSeq_sequence').text:
                source_seq_tmp.remove(seq)
                source_seq_tmp.add(sequence.find('INSDSeq_sequence').text)
                break
            else:
                source_seq_tmp.add(sequence.find('INSDSeq_sequence').text)
    source_seq = source_seq_tmp.copy()
    for sequence in ROOT.iter('INSDSeq'):
        if sequence.find('INSDSeq_sequence') is None:
            nr_no_seq += 1
            break
        for seq in source_seq:
            if sequence.find('INSDSeq_sequence').text in seq and \
                    len(sequence.find('INSDSeq_sequence').text) < len(seq):
                ROOT.remove(sequence)
                nr_degseq_rm += 1
                break
    TREE.write(ARGS.o)
    print('{} sequences not found for INSDSeq entry: '.format(nr_no_seq))
    print('{} sequences removed due to sequence degeneracy'.
          format(nr_degseq_rm))
    print('{} sequences removed due to nucleotide degeneracy'.
          format(nr_degnuc_rm))
    os.system('spd-say -l sv "Programmet klart"')
