#!/usr/bin/python3
"""Demultiplex NGS data with multiplexed PCR amplicons listed in a fasta-file
"""
from fasta import FastaList


def reduce_fa(fl):
    """removes duplicate reads and adds '_nr_of _reads to the id of the
    reduced FastaList"""
    import copy
    reddic = dict()
    reddic_new = dict()
    dnaseqs = []
    for item in fl.seq_list:
        if len(item) >= ARGS.m:
            dnaseqs.append(item.split('\n')[1])
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
    reduced_sorted = sorted(reddic_new.items(), key=lambda t: t[1], reverse=True)
    nr_red_seqs = 0
    for item in reduced_sorted:
        nr_red_seqs += item[1]
    fi = open('result.txt', 'w')
    for item in reduced_sorted:
        if item[1]/reduced_sorted[0][1] < ARGS.f:
            break
        fi.write('seq: {}, count: {}\n'. format(item[0], item[1]))
    exit(0)
    return fl


if __name__ == "__main__":
    import argparse
    import shutil
    import datetime
    import getpass
    import os

    PARSER = argparse.ArgumentParser(description='Demultiplex amplicon sequen-'
                                                 'cing NGS-data (-s) based on '
                                                 'tags listed in a separate'
                                                 'fasta-file (-t'
                                                 'Assumes the usual file '
                                                 'extensions')
    PARSER.add_argument('-sd', type=str, help='Directory containing NGS-data'
                                              'fasta/q files',
                        required=True)
    PARSER.add_argument('-t', type=str, help='tag sequence-list fasta filename',
                        required=True)
    PARSER.add_argument('-m', type=int, help='minimum sequence length',
                        default=0, required=False)
    PARSER.add_argument('-f', type=float, help='minimum fraction of total nr of'
                                               'seqs', default=0,
                        required=False)
    ARGS = PARSER.parse_args()
    tag_fa = FastaList(ARGS.t)
    logfile = open('demultiplex.log', 'w')

    logfile.write('Log for ampli_deplx.py at {}\nUser: {}\n'.format(
        str(datetime.datetime.now()).split('.')[0], getpass.getuser()))
    logfile.write('Minimum sequence length = {}\n'.format(ARGS.m))
    for seqfile in os.listdir(ARGS.sd):
        seq_fa = reduce_fa(FastaList(ARGS.sd + seqfile))
        # if os.path.isdir(ARGS.sd + seqfile.split('.')[0]):
        #     shutil.rmtree(ARGS.sd + seqfile.split('.')[0])
        os.mkdir(ARGS.sd + seqfile.split('.')[0])
        os.chdir(ARGS.sd + seqfile.split('.')[0])
        nr_seq_demult = 0
        taglist = tag_fa.seq_list
        taglist_rc = tag_fa.seq_list_revc()
        for seq in range(seq_fa.nr_seq):
            rc = True
            for tag in range(tag_fa.nr_seq):
                if not rc:
                    test_tag = taglist[tag].split('\n')[1]
                else:
                    test_tag = taglist_rc[tag].split('\n')[1]
                rc = not rc
                if test_tag in seq_fa.seq_list[seq].split('\n')[1]:
                    with open(tag_fa.id_list[tag] + '.fa', 'a') as fi:
                        if len(seq_fa.seq_list[seq]) >= ARGS.m:
                            fi.write('>' + seq_fa.id_list[seq] +
                                     seq_fa.seq_list[seq])
                            nr_seq_demult += 1
        logfile.write(("{}: {} % of {} sequences demultiplexed\n".format(
            os.path.basename(ARGS.sd + seqfile),
            round(100*nr_seq_demult/seq_fa.nr_seq), seq_fa.nr_seq)))
