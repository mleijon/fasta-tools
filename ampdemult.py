#!/usr/bin/python3
"""Demultiplex NGS data with multiplexed PCR amplicons listed in a fasta-file
"""
from fasta import FastaList


def sep(opsys):
    if opsys == 'nt':
        return '\\'
    else:
        return '/'


def reduce_fa(fl):
    """removes duplicate reads and adds '_nr_of _reads to the id of the
    reduced FastaList"""
    def rc(sequence):
        sequence = sequence[::-1]
        rc_sequence = ''
        for char in sequence:
            if char.capitalize() == 'A':
                rc_sequence += 'T'
            elif char.capitalize() == 'T':
                rc_sequence += 'A'
            elif char.capitalize() == 'G':
                rc_sequence += 'C'
            elif char.capitalize() == 'C':
                rc_sequence += 'G'
            else:
                sys.exit('Non-standard charachter in sequence. Exits')
        return rc_sequence
    import copy
    reddic_new = dict()
    dnaseqs = []
    for item in fl.seq_list:
        if len(item) >= ARGS.m:
            tlst = tag_fa.seq_list
            tlst_rc = tag_fa.seq_list_revc()
            for tagitem in range(0, tag_fa.nr_seq, 2):
                testt1 = [tlst[tagitem + 1].split('\n')[1],
                          tlst_rc[tagitem].split('\n')[1]]
                testt2 = [tlst[tagitem].split('\n')[1],
                          tlst_rc[tagitem + 1].split('\n')[1]]
                if all(x in item for x in testt1):
                    dnaseqs.append(item.split('\n')[1])
                    break
                elif all(x in item for x in testt2):
                    dnaseqs.append(rc(item.split('\n')[1]))
                    break
    for item in dnaseqs:
        reddic = copy.deepcopy(reddic_new)
        if len(reddic) == 0:
            reddic_new[item] = 1
            reddic[item] = 1
        else:
            count = 0
            for key in reddic:
                count += 1
                if item in key:
                    reddic_new[key] += 1
                    break
                elif key in item:
                    reddic_new[item] = reddic_new[key] + 1
                    del reddic_new[key]
                    break
                elif count == len(reddic):
                    reddic_new[item] = 1
    reduced_sorted = sorted(reddic_new.items(), key=lambda t: t[1],
                            reverse=True)
    fl.seq_list = []
    fl.id_list = []
    count = 0
    for item in reduced_sorted:
        if item[1] < ARGS.c:
            break
        count += 1
        fl.id_list.append('{}_{:08d}_count:{}_length:{}'.format(fl.name[fl.name.rfind(
            'Cp'):fl.name.rfind('.fa')], count, item[1], len(item[0])))
        fl.seq_list.append(('{}\n{}\n'.format(fl.id_list[count - 1], item[0])))
    fl.nr_seq = len(fl.seq_list)
    return fl


if __name__ == "__main__":
    import argparse, shutil, datetime, getpass, os, sys

    PARSER = argparse.ArgumentParser(description='Demultiplex amplicon sequen-'
                                                 'cing NGS-data (-s) based on '
                                                 'tags listed in a separate'
                                                 'fasta-file (-t'
                                                 'Assumes the usual file '
                                                 'extensions')
    PARSER.add_argument('-id', type=str, help='Directory for input data',
                        required=True)
    PARSER.add_argument('-od', type=str, help='Directory for output data',
                        required=True)
    PARSER.add_argument('-t', type=str, help='tag sequence-list fasta filename',
                        required=True)
    PARSER.add_argument('-m', type=int, help='minimum sequence length',
                        default=0, required=False)
    PARSER.add_argument('-c', type=int, help='minimum nr of seqs', default=0,
                        required=False)
    PARSER.add_argument('-f', type=float, help='minimum nr of seqs', default=0,
                        required=False)
    ARGS = PARSER.parse_args()
    if not os.path.isfile(ARGS.t):
        sys.exit('No PCR tag file. Exits.')
    if not ((ARGS.id.endswith(sep(os.name)) and
             ARGS.od.endswith(sep(os.name)))):
        sys.exit('Invalid directory name. Exits.')
    if os.path.isfile(ARGS.od[:-1]):
        print('{} is a file'.format(ARGS.od[:-1]))
        sys.exit('Exits')
    if os.path.isdir(ARGS.od):
        shutil.rmtree(ARGS.od)
    os.mkdir(ARGS.od)
    tag_fa = FastaList(ARGS.t)
    logfile = open(ARGS.od + 'demultiplex.log', 'w')
    logfile.write('Log for ampdemult.py at {}\nUser: {}\n'.format(
        str(datetime.datetime.now()).split('.')[0], getpass.getuser()))
    logfile.write('Minimum sequence length = {}\n'.format(ARGS.m))
    logfile.write('Minimum nr of sequences = {}\n'.format(ARGS.c))
    nr_of_files = len([name for name in os.listdir(ARGS.id) if
                       os.path.isfile(ARGS.id + name)])
    file_nr = 1
    for seqfile in os.listdir(ARGS.id):
        print('\rprocessing file {}/{}'.format(file_nr, nr_of_files), end=" ")
        seq_fa = reduce_fa(FastaList(ARGS.id + seqfile))
        nr_seq_demult = 0
        taglist = tag_fa.seq_list
        taglist_rc = tag_fa.seq_list_revc()
        tag_maxcount = dict()
        for item in tag_fa.id_list[::2]:
            tag_maxcount[item[:-3]] = 0
        for seq in range(seq_fa.nr_seq):
            for tag in range(0, tag_fa.nr_seq, 2):
                test_tags1 = [taglist[tag+1].split('\n')[1],
                              taglist_rc[tag].split('\n')[1]]
                test_tags2 = [taglist[tag].split('\n')[1],
                              taglist_rc[tag + 1].split('\n')[1]]
                if all(x in seq_fa.seq_list[seq] for x in test_tags1) or\
                        all(x in seq_fa.seq_list[seq] for x in test_tags2):
                    # if tag-count-max = 0 set tag-count-max = seq_fa
                    with open(ARGS.od + tag_fa.id_list[tag][:-3] + '.fa', 'a')\
                            as fi:
                        # if ARGS.f > 0
                        #seqname = seq_fa.seq_list[seq].split(\n)[0]
                        #seqcount = seqname.split(':')[1].split('_')[0]
                        #seqseq = seq_fa.seq_list[seq].split(\n)[1]
                        #newseq = '{}+_fraction:{}\n{}\n'.format(seqname, seqcount, seqseq)
                        #else:
                        fi.write('>' + seq_fa.seq_list[seq])
                        nr_seq_demult += 1
        if seq_fa.nr_seq == 0:
            seq_fa.nr_seq = 1
        logfile.write(("{}: {} % of {} sequences demultiplexed\n".format(
            os.path.basename(ARGS.id + seqfile),
            round(100 * nr_seq_demult / seq_fa.nr_seq), seq_fa.nr_seq)))
        file_nr += 1

