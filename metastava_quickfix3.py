#!/usr/bin/python
"""Metastava quickfix2 - will create a two column file with sequence name and species assignment"""

import linecache as lc
import argparse
import collections
import blast

PARSER = argparse.ArgumentParser(description='TBD')
PARSER.add_argument('-b', type=str, help='Blast hit summary table', default='test.blast')
PARSER.add_argument('-n', type=str, help='taxid name file', default='names.dmp')
PARSER.add_argument('-a', type=str, help='acc. to taxid file', default='acc2taxid.dmp')
PARSER.add_argument('-o', type=str, help='Output file', default='./nonsal.txt')
ARGS = PARSER.parse_args()
bl_inf = open(ARGS.b)
outf = open(ARGS.o, 'w')
result = dict()
acc2taxid = dict()
taxid2name = dict()


# def process_work():
#     for count, line in enumerate(bl_inf):
#         if line[0:6] == 'Query=':
#             seq_name = line.split(' ')[1].strip()
#             if 'No hits found' in lc.getline(ARGS.b, count + 7):
#                 result[seq_name] = 'NOT FOUND'
#             else:
#                 acc = lc.getline(ARGS.b, count+9).split(' ')[0].strip()[4:]
#                 acc_pointer = mtaxidf.find(acc.encode('UTF-8'))
#                 if acc_pointer == -1:
#                     result[seq_name] = 'NOT FOUND'
#                 else:
#                     mtaxidf.seek(acc_pointer)
#                     tax_id = mtaxidf.readline().decode('UTF-8').split('\t')[1].strip()
#                     mnamesf.seek(mnamesf.find(tax_id.encode('UTF-8')))
#                     result[seq_name] = mnamesf.readline().decode('UTF-8').split('|')[1].strip()
#         mnamesf.seek(0)
#         mtaxidf.seek(0)
#     return result


if __name__ == '__main__':
    with open(ARGS.a) as fi:
        for line in fi:
            acc2taxid[line.split()[1].strip()] = line.split()[2].strip()
    with open(ARGS.n) as fi:
        for line in fi:
            if line.split('|')[3].strip() == 'scientific name':
                taxid2name[line.split('|')[0].strip()] = line.split('|')[1].strip()
    print(taxid2name[0])
    print(taxid2name[1])
    print('**')
    print(taxid2name[0])
    print(taxid2name[1])
    exit(0)
    ordered_result = collections.OrderedDict(sorted(result.items()))
    for seq, spec in ordered_result.items():
        outf.write('%s\t%s\n' % (seq, spec))
    namesf.close()
    taxidf.close()
    bl_inf.close()
    outf.close()



