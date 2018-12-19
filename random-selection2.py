#!/usr/bin/python
""""Random selection among elements in an input file"""

import argparse
import random

PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-f', type=str, help='Input list of items', required=True)
PARSER.add_argument('-s', type=int, help='Sample size', default=5)
ARGS = PARSER.parse_args()
itemlist = []
samplelist = []
with open(ARGS.f) as fi:
    for line in fi:
        itemlist.append(tuple(line.split()))
while len(samplelist) < ARGS.s:
    if len(itemlist) == 0:
        print('Can only sample %d items' % len(samplelist))
        exit(0)
    else:
        sample = random.sample(itemlist,1)[0]
        samplelist.append(sample)
        itemlist[:] = [x for x in itemlist if not x[1] == sample[1]]
print(samplelist)
###