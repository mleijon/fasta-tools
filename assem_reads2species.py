import argparse
from collections import OrderedDict
PARSER = argparse.ArgumentParser(description='Aseemble reads to species.')
PARSER.add_argument('-f', type=str, help='read list file', required=True)
PARSER.add_argument('-o', type=str, help='read list file', required=True)
ARGS = PARSER.parse_args()
with open(ARGS.f) as fi:
    inpdata = fi.read()
    readlst = []
    for item in inpdata.split('\n')[:-1]:
        if ' ' in item:
            readlst.append(item.split('\t')[1].split(' ')[0] + ' '
                           + item.split('\t')[1].split(' ')[1])
        else:
            readlst.append(item.split('\t')[1])
readset = set(readlst)
readdict = {key: round(100*inpdata.count(key)/len(readlst), 1) for key in
            readset}
orderedreads = OrderedDict(sorted(readdict.items(), key=lambda t: t[1],
                                  reverse=True))
with open(ARGS.o, 'w') as fi:
    for key in orderedreads:
        if orderedreads[key] > 1:
            fi.write(key + '\t\t' + str(orderedreads[key]) + ' %\n')
    fi.write('\nTotal number of reads: {}'.format(len(readlst)))


