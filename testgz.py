#!/usr/bin/python3
import gzip
outfile = open('U.fastq','w')
with gzip.open('U.fastq.gz','rt') as f:
    test = f.read()
outfile.write(test)
outfile.close()
