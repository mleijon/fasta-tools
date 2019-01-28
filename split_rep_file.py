#!/usr/bin/python
"""Split a repeating file in n smaller files"""

import tempfile as tf


def split(infi, sep, nr_of_splits):
    tmp_file_list = []
    for count in range(nr_of_splits):
        tmp_file_list.append(tf.TemporaryFile('w+t'))
    nr_of_elem = 0
    for line in infi:
        if sep in line:
            nr_of_elem += 1
    nr_of_elem_per_split = nr_of_elem // nr_of_splits
    nr_1_extra = nr_of_elem % nr_of_splits
    file_count = 0
    element_count = 0
    infi.seek(0)
    for line_count, line in enumerate(infi):
        if sep in line:
            elem_count += 1
        if count % (nr_of_elem_per_split + (nr_1_extra <= 0)) == 0:
            nr_1_extra -= 1
            file_count += 1






