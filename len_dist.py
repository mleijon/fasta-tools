import os
import argparse
from collections import defaultdict

PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-d', type=str, help='Input directory', required=True)
ARGS = PARSER.parse_args()
ord_counts_R1 = defaultdict(lambda: 0)
ord_counts_R2 = defaultdict(lambda: 0)
row1 = row2 = row3 = ''
all_keys = set()
all_sorted_keys = list()
seq = False
for file in [x for x in os.listdir(ARGS.d) if x.endswith('_R1_001_trimmed.fastq')
                                              and not os.path.isdir(os.path.join(ARGS.d, x))]:
    with open(os.path.join(ARGS.d, file)) as f:
        for line in f:
            if line[0] == '@':
                seq = True
                continue
            elif seq:
                ord_counts_R1[len(line)] += 1
                seq = False
        f.close()
    seq = False
    with open(os.path.join(ARGS.d, file.replace('_R1_', '_R2_'))) as f:
        for line in f:
            if line[0] == '@':
                seq = True
                continue
            elif seq:
                ord_counts_R2[len(line)] += 1
                seq = False
        f.close()
    all_keys = set(ord_counts_R1.keys())
    all_keys = all_keys.union(set(ord_counts_R2.keys()))
    all_sorted_keys = list(all_keys)
    all_sorted_keys.sort(reverse=True)
    for item in all_sorted_keys:
        row1 += str(item) + ','
        row2 += str(ord_counts_R1[item]) + ','
        row3 += str(ord_counts_R2[item]) + ','
    row1 = row1[:-1]
    row2 = row2[:-1]
    row3 = row3[:-1]
    out_file = os.path.join(ARGS.d, file.split('_')[0] + '.csv')
    with open(out_file, 'w') as f_out:
        f_out.write(row1 + '\n')
        f_out.write(row2 + '\n')
        f_out.write(row3 + '\n')
    ord_counts_R1 = defaultdict(lambda: 0)
    ord_counts_R2 = defaultdict(lambda: 0)
    row1 = row2 = row3 = ''
    all_keys = set()
    all_sorted_keys = list()
    seq = False
