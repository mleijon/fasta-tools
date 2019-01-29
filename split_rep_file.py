#!/usr/bin/python
"""Split a repeating file in n smaller files"""

import tempfile as tf
tmp_file_list = []


def split(infi, sep, nr_of_splits):
    for count in range(nr_of_splits):
        tmp_file_list.append(tf.TemporaryFile('w+t'))
    nr_of_elem = 0
    for line in infi:
        if sep in line:
            nr_of_elem += 1
    nr_of_elem_per_split = nr_of_elem // nr_of_splits
    nr_1_extra = nr_of_elem % nr_of_splits
    file_count = -1
    elem_count = -1
    infi.seek(0)
    for line in infi:
        if sep in line:
            elem_count += 1
        if elem_count % (nr_of_elem_per_split + (nr_1_extra >= 0)) == 0 and\
                sep in line:
            nr_1_extra -= 1
            file_count += 1
            tmp_file_list[file_count].write(line)
        else:
            tmp_file_list[file_count].write(line)
    for count in range(nr_of_splits):
        tmp_file_list[count].seek(0)


in_fi = open('test.fa')
split(in_fi, '>', 2)
for i in range(2):
    print('i: {}'.format(i))
    print(tmp_file_list[i].read())










