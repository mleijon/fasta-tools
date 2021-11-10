#!/usr/bin/python

import split_file as sf


infile = open('test.blast')
tempfiles = sf.split(infile, "BLASTN 9.1.0", 6)
for count, item in enumerate(tempfiles):
    with open('test_' + str(count) + '.blast', 'w') as of:
        print(item.name)
        of.write(item.read())




