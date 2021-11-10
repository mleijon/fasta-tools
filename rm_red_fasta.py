from fasta import FastaList
from collections import defaultdict
import argparse
from datetime import datetime

MAX_DEG = 0
MIN_LEN = 27000
deg_counters = defaultdict(lambda: 0)


# Counts the number of deg nucs in seq
def cnt_nt_deg(seq):
    cnt = 0
    for char in seq:
        if char not in ['A', 'C', 'G', 'T']:
            cnt += 1
    return cnt


def maxmin_len(seqs):
    seqs = list(seqs)
    max_length = len(seqs[0])
    min_length = len(seqs[0])
    for seq in seqs:
        if len(seq) > max_length:
            max_length = len(seq)
        if len(seq) < min_length:
            min_length = len(seq)
    return min_length, max_length


# Creates a global dict with the nr of deg nucs as keys and the number of seqs
# carrying the number of deg nucs as value. Writes to an optional file (ARGS.d).
def deg_summary(seqs):
    global deg_counters
    for seq in seqs:
        deg_counters[cnt_nt_deg(seq)] += 1
    with open(ARGS.d, 'w') as deg_file:
        for key in sorted(deg_counters.keys()):
            deg_file.write('{}\t{}\n'.format(key, deg_counters[key]))


def red_uniq_seq(uniq_seq):
    global start
    if ARGS.subsample:
        uniq_seq = sample(uniq_seq, 10*ARGS.subsample)
        print('Subsampled: {:38}'.format(len(uniq_seq)))
    new_seq = [uniq_seq[0]]
    i = 1
    one_percent = (len(uniq_seq) + len(uniq_seq) % 100) // 100
    for s1 in uniq_seq[1:]:
        if i % one_percent == 0:
            print('\r{}%'.format(i // one_percent), end="")
        i += 1
        add = True
        for s2 in new_seq:
            if len(s1) == len(s2):
                continue
            if len(s1) < len(s2) and s1 in s2:
                add = False
                break
            elif len(s2) < len(s1) and s2 in s1:
                new_seq.remove(s2)
                new_seq.append(s1)
                add = False
                break
        if add:
            new_seq.append(s1)
    end = datetime.now()
    delta_time = end - start
    calc_time = delta_time.seconds
    print('\rDone in ', end="")
    if calc_time >= 86400:
        print('{} days'.format(calc_time // 86400), end="")
        calc_time = calc_time % 86400
    if calc_time >= 3600:
        print(' {} hours'.format(calc_time // 3600), end="")
        calc_time = calc_time % 3600
    if calc_time >= 60:
        print(' {} minutes'.format(calc_time // 60), end="")
        calc_time = calc_time % 60
    else:
        print(' {} seconds'.format(calc_time))
        print()
    return new_seq


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-i', type=str, help='input fasta file', required=True)
    PARSER.add_argument('-o', type=str, help='output fasta file', required=True)
    PARSER.add_argument('-d', type=str, help='output degeneracy file', required=False)
    PARSER.add_argument('-l', type=int, help='fasta row length', required=False)
    PARSER.add_argument('--rm_subseqs',  action='store_true', help='switch for removal of subseqs')
    PARSER.add_argument('--subsample',  type=int, help='switch for removal of subseqs', required=False)

    ARGS = PARSER.parse_args()
    if ARGS.subsample:
        from random import sample
    org_seqs = list()
    unique_seqs = set()
    #  Create a set of unique seqs without flanking Ns or polyA-tail and store
    #  the originals seqs in the list org_seqs. The seqs must be at least of
    #  MIN_LEN length and contain at most MAX_DEG number of deg nucs
    for item in FastaList(ARGS.i):
        org_seqs.append(item.split('\n')[1])
        sequence = item.split('\n')[1].strip('N').rstrip('A')
        if len(sequence) >= MIN_LEN and cnt_nt_deg(sequence) <= MAX_DEG:
            unique_seqs.add(sequence)
    if ARGS.d:
        deg_summary(org_seqs)
    print('Nr of input sequences: {:38}'.format(len(org_seqs)))
    print('Nr of unique sequences (including subsequences): {:12}'.format(len(unique_seqs)))
    unique_seqs = list(unique_seqs)
    if ARGS.rm_subseqs:
        print('\n***Removing subsequences***')
        unique_seqs.sort(key=len, reverse=True)
        start = datetime.now()
        final_seqs = red_uniq_seq(unique_seqs)
        print('Nr of unique sequences (excluding subsequenes): {:13}'.format(len(final_seqs)))
    else:
        print('subsequences not removed ')
        final_seqs = unique_seqs
    count = 0
    if ARGS.subsample:
        final_seqs = sample(final_seqs, ARGS.subsample)
        print('Subsampled: {:38}'.format(len(final_seqs)))
    else:
        print('No subsampling')
    with open(ARGS.o, 'w') as fi:
        for s in final_seqs:
            count += 1
            if ARGS.l:
                s_list = [s[i:i + ARGS.l] for i in range(0, len(s), ARGS.l)]
                s = '\n'.join(s_list)
            fi.write('>Uniq_{}\n{}\n'.format(count, s))
