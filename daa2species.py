#!/usr/bin/python3

sciNames = dict()
mergedTaxids = dict()
deletedTaxids = set()
names_file = open("names.dmp")
merge_file = open("merged.dmp")
delete_file = open("delnodes.dmp")
for line in merge_file:
    item = line.split('|')
    mergedTaxids.update({item[0].strip(): item[1].strip()})
for line in delete_file:
    item = line.split('|')[0].strip()
    deletedTaxids.add(item)
merge_file.close()
daa_file = open("SE-2194-Botvid-Granby-filter.daa")
for line in names_file:
    item = line.split('|')
    if "scientific" in item[3]:
        sciNames.update({item[0].strip(): item[1].strip()})
names_file.close()
for line in daa_file:
    item = line.split()[1].strip()
    if item != '0':
        if item in mergedTaxids.keys():
            print(sciNames[mergedTaxids[item]])
        elif item in deletedTaxids:
            print(item + ' :deleted')
        else:
            print(sciNames[item])
            
daa_file.close()
