#!/usr/bin/python
FASTQ_EXT = ".fastq"
import os
import argparse

PARSER = argparse.ArgumentParser(description='Calculate average read length for all fastq-files in a directory')
PARSER.add_argument('-d', type=str, help='file directory file', required=True)
ARGS = PARSER.parse_args()
nr_of_seq = 0
length_sum = 0
file_count = 0
for filename in os.listdir(ARGS.d):
    if filename.endswith(FASTQ_EXT):
        prev_is_lbl = False
        file_count +=1
        with open(os.path.join(ARGS.d, filename)) as f:
            print('Analysing file: {}'.format(filename))
            for line in f:
                if prev_is_lbl:
                    length_sum += len(line)
                    nr_of_seq += 1
                    prev_is_lbl = False
                elif line[0] == '@':
                    prev_is_lbl = True
print(nr_of_seq)
print('Average read length is {}'.format(round(length_sum/nr_of_seq), 0))
print('Average nr of reads is {}'.format(round(3*nr_of_seq/file_count), 0))
