#!/usr/bin/python
"""This """

import json

readfi = open('foodPT2018.fastq')
salfi = open('salhits.txt')
nosalfi = open('non_salhits.txt')
svaresult = open('svaresult.txt', 'w')
namelist = []
sallist = []
nosallist = []

try:
    with open('foodPT2018.js') as name_file:
        namelist = json.load(name_file)
except FileNotFoundError:
    for line in readfi:
        if line[0:8] == '@fpt2018':
            namelist.append(line[8:].lstrip('0').strip())
    with open('foodPT2028.js', 'w') as cs_file:
        json.dump(namelist, cs_file, indent=4)
for line in salfi:
    sallist.append(line.split()[0][7:].lstrip('0'))
for line in nosalfi:
    nosallist.append(line.split()[0][7:].lstrip('0'))
salfi.seek(0)
nosalfi.seek(0)
current_nosal = nosalfi.readline()
current_nosal_nr = int(current_nosal.split()[0][7:].lstrip('0'))
current_sal = salfi.readline()
current_sal_nr = int(current_sal.split()[0][7:].lstrip('0'))
eof_sal = False
eof_nosal = False
for count in namelist:
    if count == sallist[0] and not eof_sal:
        svaresult.write(current_sal)
        _ = sallist.pop(0)
        current_sal = salfi.readline()
        if current_sal == '':
            eof_sal = True
        else:
            current_sal_nr = int(current_sal.split()[0][7:].lstrip('0'))
        continue
    if count == nosallist[0] and not eof_nosal:
        svaresult.write(current_nosal)
        _ = nosallist.pop(0)
        current_nosal = nosalfi.readline()
        if current_nosal == '':
            eof_nosal = True
        else:
            current_nosal_nr = int(current_nosal.split()[0][7:].lstrip('0'))
        continue
    if (count != sallist[0]) and (count != nosallist[0]):
        svaresult.write('fpt2018'+str(count).zfill(8)+'\tNo hit\n')
    print(len(sallist))
    print(len(nosallist))
salfi.close()
nosalfi.close()
svaresult.close()
