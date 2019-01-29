#!/usr/bin/python
"""Split a repeating file in n smaller files. The file must contain repeating
occurences of a tag sequence, for instance the '>' of a fasta-file. In addition
the first line of the file must contain the tag. The splitted files will all
start with a line containing the tag"""

import tempfile as tf
tmp_file_list = []


def split(infi, sep, nr_of_splits):
    nr_of_elem = 0
    for line in infi:
        if sep in line:
            nr_of_elem += 1
    nr_of_elem_per_split = nr_of_elem // nr_of_splits
    if nr_of_elem_per_split == 0:  # Handles nr_of_splits > nr_pf_elem
        nr_of_elem_per_split = 1
        nr_of_splits = nr_of_elem
        nr_1_extra = 0
    else:
        nr_1_extra = nr_of_elem % nr_of_splits
    for count in range(nr_of_splits):
        tmp_file_list.append(tf.TemporaryFile('w+t'))
    file_count = 0
    elem_count = 0
    infi.seek(0)
    for line in infi:
        if sep in line:
            elem_count += 1
        if elem_count % (nr_of_elem_per_split + (nr_1_extra > 0)) == 0 and\
                sep in line:
            print(nr_of_elem_per_split + (nr_1_extra > 0))
            nr_1_extra -= 1
            tmp_file_list[file_count].write(line)
            file_count += 1
        else:
            tmp_file_list[file_count].write(line)
    for count in range(nr_of_splits):
        tmp_file_list[count].seek(0)


in_fi = open('DNA-VIR_S2.fa')
split(in_fi, '>', 12)
for i in range(len(tmp_file_list)):
    filename = 'DNA-VIR_S2' + '_' + str(i) + '.fa'
    fi = open(filename, 'w')
    fi.write(tmp_file_list[i].read())
    fi.close()











