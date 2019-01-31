#!/usr/bin/python
"""Split a repeating file in n smaller temporary files.

The parameters of the function is the input file object 'infi', the separator
string 'sep' and the number of splits (files) 'n' created. The file must contain
repeating occurences of a separator sequence, for instance ,the '>' of a
fasta-file. The splitted files will all start with a line containing the
separator, except the first file that may start with lines without the
separator. The lines from one separator up to but not including the next line
 with a separatoris an element (again with the first element as a possible
exception. The number of elements in the file is 'nr_of_elem'. The number of
elements per files after splitting is 'nr_of_elem_per_split'. If the nr_of_elem
are not divisible by n (number of splits) the 'nr_1_extra' first files will
contain one element more than the remaining files"""

import tempfile as tf
tmp_file_list = []


def split(infi, sep, nr_of_splits):
    nr_of_elem = 0
    for line in infi:
        if sep in line:
            nr_of_elem += 1
    if nr_of_elem == 0:
        return tmp_file_list  # Return empty list if separator not present
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
        if elem_count > nr_of_elem_per_split + (nr_1_extra > 0):
            elem_count = 1
            file_count += 1
            nr_1_extra -= 1
        tmp_file_list[file_count].write(line)
    for count in range(nr_of_splits):
        tmp_file_list[count].seek(0)
    return tmp_file_list  # Return a list of n temporary file objects


if __name__ == "__main__":
    import argparse as ap
    PARSER = ap.ArgumentParser(description='Splits file in smaller based\
    on a separator string')
    PARSER.add_argument('-f', type=str, help='target file', required=True)
    PARSER.add_argument('-s', type=str, help='separatr string', required=True)
    PARSER.add_argument('-n', type=int, help='nr of splits', required=True)
    ARGS = PARSER.parse_args()
    try:
        infi = open(ARGS.f)
    except FileNotFoundError:
        print('File not found. Exits.')
        exit(0)
    split(infi, ARGS.s, ARGS.n)
    for i in range(ARGS.n):
        outfi = open(infi.name.partition('.')[0] + '_' + str(i) +
                     infi.name.partition('.')[1] + infi.name.partition('.')[2], 'w')
        outfi.write(tmp_file_list[i].read())
        outfi.close()












