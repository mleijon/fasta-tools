#!/usr/bin/python

import sys

File = open(sys.argv[1])
while 1:
    Line = File.readline()
    if len(Line) == 0:
        break
    s = ""
    for c in Line[:-1]:
        if ord(c) > 127:
            print(Line[:-1])
            break
