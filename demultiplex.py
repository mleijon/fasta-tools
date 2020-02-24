#!/usr/bin/python3
"""Demultiplex NGS data (fa or fastq - may be gzipped) with primers listed in a
separate fasta-file
"""
from fasta import FastaList
import argparse
import shutil
import datetime
import getpass
import os
import sys
import copy


def sep(opsys):
    """Handles the different directory separators in win/linux
     TODO: check if necessary """
    if opsys == 'nt':
        return '\\'
    else:
        return '/'


def reduce_fa(fl):
    """removes duplicate reads and adds ' the read counts to the id of the
    reduced FastaList"""
    # TODO: Include which primerpair each sequence originate from

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

    # Loop over all sequences in the input file (in FastaList format) and check
    # that the sequence is longer than the minimum length (ARGS.m) and that a
    # primer pair is present in the sequence. Primer pairs are checked in both
    # directions. The resulting sequnces are stored without id in the list
    # dnaseqs
    dnaseqs = []  # Sequence list containing complete for and rev primers
    for item in fl.seq_list:
        if len(item.split('\n')[1]) >= ARGS.m:  # filter on seq length
            primerlst = primer_fa.seq_list
            primerlst_rc = primer_fa.seq_list_revc()
            for primernr in range(0, primer_fa.nr_seq, 2):
                # Test primers in both orientations
                testpair1 = [primerlst[primernr + 1].split('\n')[1],
                             primerlst_rc[primernr].split('\n')[1]]
                testpair2 = [primerlst[primernr].split('\n')[1],
                             primerlst_rc[primernr + 1].split('\n')[1]]
                # Ensure that both primers are in sequence
                if all(x in item for x in testpair1):
                    dnaseqs.append(item.split('\n')[1])
                    break
                elif all(x in item for x in testpair2):
                    dnaseqs.append(rc(item.split('\n')[1]))
                    break
    # Merge sequences that are identical or substrings of each other and store
    # the longest representative. Keep a count of the number of sequences that
    # are merged. Store the sequenes and the count in  sorted list of tuples:
    # (sequence, count)
    merged_seqs = dict()  # {sequence:count}
    for item in dnaseqs:
        merged_seqs_tmp = copy.deepcopy(merged_seqs)
        if len(merged_seqs_tmp) == 0:
            # Initilize count for sequence "item" when merged_seqs_tmp is empty
            merged_seqs[item] = 1
            merged_seqs_tmp[item] = 1
        else:
            count = 0
            # loop over sequences (keys) and update count
            for key in merged_seqs_tmp:
                count += 1
                if item in key:
                    merged_seqs[key] += 1
                    break
                elif key in item:
                    merged_seqs[item] = merged_seqs[key] + 1  # Keep the longest
                    del merged_seqs[key]                     # sequence
                    break
                elif count == len(merged_seqs_tmp):
                    merged_seqs[item] = 1
    # Create a sorted list of tuples (key, value) from dict merged_seqs
    # larger to smaller counts
    reduced_sorted = sorted(merged_seqs.items(), key=lambda t: t[1],
                            reverse=True)  # Sort sequences on count
    fl.seq_list = []
    fl.id_list = []
    count = 0 # Generic counter added to id string to make id:s unique
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


def main():
    # Writes log-file
    logfile = open(ARGS.od + os.path.basename(__file__).
                   replace('.py', '.log'), 'w')
    logfile.write('Log for: {} at {}\nUser: {}\n\n'.format(os.path.basename(
        __file__), str(datetime.datetime.now()).split('.')[0],
        getpass.getuser()))
    logfile.write('Minimum sequence length = {}\n'.format(ARGS.m))
    logfile.write('Minimum nr of sequences = {}\n'.format(ARGS.c))
    logfile.write('Minimum fraction of most abundant sequence = {}\n\n'.format(
        ARGS.f))

    # Loop over files in input directory (ARGS.id) but skip files without .fa
    # and .fastq file extension
    filelst = [name for name in os.listdir(ARGS.id) if
               os.path.isfile(ARGS.id + name) and (name.endswith('.fa') or
                                                   name.endswith('.fastq'))]
    nr_of_files = len(filelst)
    file_nr = 1
    for seqfile in filelst:
        print('\rprocessing file {}/{}'.format(file_nr, nr_of_files), end=" ")
        inp_seq = FastaList(ARGS.id + seqfile)
        init_seq = inp_seq.nr_seq  # The intitial number of seqs in fasta-file
        seq_fa = reduce_fa(inp_seq)
        logfile.write('{}: Read {} sequences. '.format(seqfile, init_seq))
        nr_seq_demult = 0
        primerlist = primer_fa.seq_list
        primerlist_rc = primer_fa.seq_list_revc()
        marker_maxcount = dict()
        fraction = 0
        if ARGS.f > 0:
            # Assumes three charachters at the end of the primer id indicating
            # forward or reverse. Initilize the dict contaning the counts for
            # the sequence for each primer pair with the  highest count
            # TODO: Try to remove dependence on specific primer names in
            #  primer-file
            for primerid in primer_fa.id_list[::2]:
                marker_maxcount[primerid[:-3]] = 0
                
        for seq in range(seq_fa.nr_seq):
            seqname = seq_fa.seq_list[seq].split('\n')[0]
            seqcount = int(seqname.split(':')[1].split('_')[0])
            seqseq = seq_fa.seq_list[seq].split('\n')[1]
            for primer in range(0, primer_fa.nr_seq, 2):
                test_primers1 = [primerlist[primer+1].split('\n')[1],
                                 primerlist_rc[primer].split('\n')[1]]
                test_primers2 = [primerlist[primer].split('\n')[1],
                                 primerlist_rc[primer + 1].split('\n')[1]]
                if all(x in seq_fa.seq_list[seq] for x in test_primers1) or\
                        all(x in seq_fa.seq_list[seq] for x in test_primers2):
                    if ARGS.f > 0:
                        # Since seq_fa.seq_list is ordered with respect to count
                        # the first occurance in the list of a marker sequence
                        # will have the highest count
                        if marker_maxcount[primerlist[primer].split('\n')[0][1:-3]] == 0:
                            marker_maxcount[primerlist[primer].split('\n')[0][1:-3]] =\
                                seqcount
                        fraction = seqcount / marker_maxcount[
                            primerlist[primer].split('\n')[0][1:-3]]
                        # Filters on fraction of most abundant sequnce for the
                        # marker
                        if fraction < ARGS.f:
                            break
                    with open(ARGS.od + primer_fa.id_list[primer][:-3] + '.fa', 'a')\
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


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='The script demultiplex '
                                                 'amplicon sequencing NGS-data'
                                                 ' based on primers (multiplexed '
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
    PARSER.add_argument('-t', type=str, help='fasta file name for primer/primer '
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
        sys.exit('No PCR primer file. Exits.')
    if not ((ARGS.id.endswith(sep(os.name)) and
             ARGS.od.endswith(sep(os.name)))):
        sys.exit('Invalid directory name. Exits.')
    if os.path.isfile(ARGS.od[:-1]):
        print('{} is a file'.format(ARGS.od[:-1]))
        sys.exit('Exits')
    if os.path.isdir(ARGS.od):
        shutil.rmtree(ARGS.od)
    os.mkdir(ARGS.od)
    primer_fa = FastaList(ARGS.t)  # Make fastaList of primerfile
    if not (0 <= ARGS.f <= 1):
        sys.exit('Fraction (-f) out of range. Exits.')
    if ARGS.m < 0:
        sys.exit('Minimum sequence length (-m) ut of range. Exits.')
    if ARGS.c < 1:
        sys.exit('Count (-c) ot of range. Exits.')
    main()
