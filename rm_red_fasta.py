from fasta import FastaList
from collections import defaultdict
import argparse
from datetime import datetime

MAX_DEG = 0
MIN_LEN = 27000

PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-i', type=str, help='input fasta file', required=True)
PARSER.add_argument('-o', type=str, help='output fasta file', required=True)
PARSER.add_argument('-d', type=str, help='output degeneracy file', required=False)
PARSER.add_argument('-l', type=int, help='fasta row length', required=False)
ARGS = PARSER.parse_args()
unique_seqs = set()
org_seqs = list()
deg_counters = defaultdict(lambda: 0)
fl = FastaList(ARGS.i)
max_len = 0
min_len = 27000


def cnt_nt_deg(seq):
    cnt = 0
    for char in seq:
        if char not in ['A', 'C', 'G', 'T']:
            cnt += 1
    return cnt


for item in fl:
    org_seqs.append(item.split('\n')[1])
    sequence = item.split('\n')[1].strip('N').rstrip('A')
    if len(sequence) >= MIN_LEN and cnt_nt_deg(sequence) <= MAX_DEG:
        unique_seqs.add(sequence)


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


def deg_summary(seqs):
    global deg_counters
    for seq in seqs:
        deg_counters[cnt_nt_deg(seq)] += 1
    if ARGS.d:
        with open(ARGS.d, 'w') as deg_file:
            for key in sorted(deg_counters.keys()):
                deg_file.write('{}\t{}\n'.format(key, deg_counters[key]))


def red_uniq_seq(uniq_seq):
    global start
    new_seq = [uniq_seq[0]]
    i = 1
    one_percent = (len(uniq_seq) + len(uniq_seq) % 100) // 100
    for s1 in uniq_seq:
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
        print(' {} seconds'.format(calc_time))
        print()
    return new_seq


deg_summary(org_seqs)
print('Nr of input sequences: {:38} {}'.format(len(org_seqs), maxmin_len(org_seqs)))
print('Nr of sequences without degeneracy: {:25}'.format(deg_counters[0]))
print('Nr of unique sequences (including subsequences): {:12} {}'.format(len(unique_seqs), maxmin_len(unique_seqs)))
print('\n***Removing subsequences***')
unique_lst = list(unique_seqs)
unique_lst.sort(key=len, reverse=True)
start = datetime.now()
final_seqs = red_uniq_seq(unique_lst)
print('Nr of unique sequences (excluding subsequenes): {:13} {}'.format(len(final_seqs), maxmin_len(final_seqs)))
count = 0
with open(ARGS.o, 'w') as fi:
    for s in final_seqs:
        count += 1
        if ARGS.l:
            s_list = [s[i:i + ARGS.l] for i in range(0, len(s), ARGS.l)]
            s = '\n'.join(s_list)
        fi.write('>Uniq_{}\n{}\n'.format(count, s))
