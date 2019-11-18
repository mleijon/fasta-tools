#!/usr/bin/python3
"""Demultiplex NGS data (fa or fastq - may be gzipped) with primers listed in a
separate fasta-file
"""
from fasta import FastaList


def sep(opsys):
    """Handles the different directory separators in win/linux
     TODO: check if necessary """
    if opsys == 'nt':
        return '\\'
    else:
        return '/'


def reduce_fa(fl):
    """removes duplicate reads and adds '_nr_of _reads to the id of the
    reduced FastaList"""
    def rc(sequence):
        """Determines and returns the reverse complement of an input sequence
        Only A,C,G and T allowed"""
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
    dnaseqs = []  # Lists sequences containing both complete for and rev primers
    for item in fl.seq_list:
        if len(item) >= ARGS.m:  # filter on seq length
            tlst = tag_fa.seq_list
            tlst_rc = tag_fa.seq_list_revc()
            for tagitem in range(0, tag_fa.nr_seq, 2):
                # Test primers in both orientations
                testt1 = [tlst[tagitem + 1].split('\n')[1],
                          tlst_rc[tagitem].split('\n')[1]]
                testt2 = [tlst[tagitem].split('\n')[1],
                          tlst_rc[tagitem + 1].split('\n')[1]]
                # Ensure that both primers are in sequence
                if all(x in item for x in testt1):
                    dnaseqs.append(item.split('\n')[1])
                    break
                elif all(x in item for x in testt2):
                    dnaseqs.append(rc(item.split('\n')[1]))
                    break
    for item in dnaseqs:
        reddic = copy.deepcopy(reddic_new)  # Carries count for sequences
        if len(reddic) == 0:
            reddic_new[item] = 1  # Initilize count for sequence "item"
            reddic[item] = 1      # when reddic is empty
        else:
            count = 0
            for key in reddic:  # loop over sequences (keys) and update count
                count += 1
                if item in key:
                    reddic_new[key] += 1
                    break
                elif key in item:
                    reddic_new[item] = reddic_new[key] + 1  # Keep the longest sequence
                    del reddic_new[key]
                    break
                elif count == len(reddic):
                    reddic_new[item] = 1
    # Create a list of tuples (key, value) from dict reddic_new
    reduced_sorted = sorted(reddic_new.items(), key=lambda t: t[1],
                            reverse=True)  # Sort sequences on count
    fl.seq_list = []
    fl.id_list = []
    count = 0
    for item in reduced_sorted:
        if item[1] < ARGS.c:  # filter sequences on count
            break
        count += 1
        # TODO: Remove assumption on filename starting with "Cp"
        fl.id_list.append('{}_{:08d}_count:{}_length:{}'.format(fl.name[fl.name.find(
            'Cp'):fl.name.rfind('.fa')], count, item[1], len(item[0])))
        fl.seq_list.append(('{}\n{}\n'.format(fl.id_list[count - 1], item[0])))
    fl.nr_seq = len(fl.seq_list)
    return fl


if __name__ == "__main__":
    import argparse
    import shutil
    import datetime
    import getpass
    import os
    import sys

    PARSER = argparse.ArgumentParser(description='The script demultiplex '
                                                 'amplicon sequencing NGS-data'
                                                 ' based on tags (multiplexed '
                                                 'primers) listed in a separate'
                                                 ' fasta-file (alternating '
                                                 'forward and reverse primers).'
                                                 ' Assumes the usual file '
                                                 'extensions. Fasta- and fastq-'
                                                 'files that may be gzipped are'
                                                 ' allowed. Sequence reads must'
                                                 ' contain complete forward '
                                                 'and reverse primer and are '
                                                 'optinally filtered on length,'
                                                 ' count (multipicity) and '
                                                 'fraction af the most '
                                                 'abundant sequence.')
    PARSER.add_argument('-id', type=str, help='Directory for input data',
                        required=True)
    PARSER.add_argument('-od', type=str, help='Directory for output data',
                        required=True)
    PARSER.add_argument('-t', type=str, help='fasta file name for tag/primer '
                                             'sequence-list',
                        required=True)
    PARSER.add_argument('-m', type=int, help='minimum seq length',
                        default=250, required=False)
    PARSER.add_argument('-c', type=int, help='minimum nr of seqs', default=1,
                        required=False)
    PARSER.add_argument('-f', type=float, help='minimum fraction of most '
                                               'abundant seq', default=0,
                        required=False)
    ARGS = PARSER.parse_args()
    # Some control of input file/directory names and parameter values
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
    if not (0 <= ARGS.f <= 1):
        sys.exit('Fraction (-f) out of range. Exits.')
    if ARGS.m < 0:
        sys.exit('Minimum sequence length (-m) ut of range. Exits.')
    if ARGS.c < 1:
        sys.exit('Count (-c) ot of range. Exits.')
    # Writes log-file
    logfile = open(ARGS.od + 'demultiplex.log', 'w')
    logfile.write('Log for ampdemult.py at {}\nUser: {}\n\n'.format(
        str(datetime.datetime.now()).split('.')[0], getpass.getuser()))
    logfile.write('Minimum sequence length = {}\n'.format(ARGS.m))
    logfile.write('Minimum nr of sequences = {}\n'.format(ARGS.c))
    logfile.write('Minimum fraction of most abundant sequence = {}\n\n'.format(
        ARGS.f))
    # Loop over files in input directory (ARGS.id)
    nr_of_files = len([name for name in os.listdir(ARGS.id) if
                       os.path.isfile(ARGS.id + name)])
    file_nr = 1
    for seqfile in os.listdir(ARGS.id):
        print('\rprocessing file {}/{}'.format(file_nr, nr_of_files), end=" ")
        inp_seq = FastaList(ARGS.id + seqfile)
        init_seq = inp_seq.nr_seq
        seq_fa = reduce_fa(inp_seq)
        logfile.write('{}: Read {} sequences. '.format(
            os.path.basename(ARGS.id + seqfile), init_seq))
        nr_seq_demult = 0
        taglist = tag_fa.seq_list
        taglist_rc = tag_fa.seq_list_revc()
        tag_maxcount = dict()
        fraction = 0
        if ARGS.f > 0:
            # TODO: Try to remove dependence on specific primer names in
            #  tag-file
            for tagid in tag_fa.id_list[::2]:  # Tag-file contain for+rev primer
                tag_maxcount[tagid[:-3]] = 0   # Assumes three charachters at
        for seq in range(seq_fa.nr_seq):       # the end indicating for/rev
            seqname = seq_fa.seq_list[seq].split('\n')[0]
            seqcount = int(seqname.split(':')[1].split('_')[0])
            seqseq = seq_fa.seq_list[seq].split('\n')[1]
            for tag in range(0, tag_fa.nr_seq, 2):
                test_tags1 = [taglist[tag+1].split('\n')[1],
                              taglist_rc[tag].split('\n')[1]]
                test_tags2 = [taglist[tag].split('\n')[1],
                              taglist_rc[tag + 1].split('\n')[1]]
                if all(x in seq_fa.seq_list[seq] for x in test_tags1) or\
                        all(x in seq_fa.seq_list[seq] for x in test_tags2):
                    if ARGS.f > 0:
                        if tag_maxcount[taglist[tag].split('\n')[0][1:-3]] == 0:
                            tag_maxcount[taglist[tag].split('\n')[0][1:-3]] =\
                                seqcount
                        fraction = seqcount / tag_maxcount[
                            taglist[tag].split('\n')[0][1:-3]]
                        if fraction < ARGS.f:  # Filters on fraction of most
                            break              # abundant
                    with open(ARGS.od + tag_fa.id_list[tag][:-3] + '.fa', 'a')\
                            as fi:
                        if ARGS.f > 0:
                            newseq = '>{}_fraction:{}\n{}\n'.format(
                                seqname, round(fraction, 3), seqseq)
                            fi.write(newseq)
                            nr_seq_demult += 1
                        else:
                            fi.write('>' + seq_fa.seq_list[seq])
                            nr_seq_demult += 1
                    fi.close()
        if seq_fa.nr_seq == 0:
            seq_fa.nr_seq = 1
        logfile.write("Reduced to {} sequences.\n".format(nr_seq_demult))
        file_nr += 1

