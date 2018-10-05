#!/usr/bin/python3
#counts the rows in infile
def countlines(infile):
    with infile as f:
        for count, element in enumerate(f,1):
            pass
    return count
