#!/usr/bin/python
species_list = []
nor_list = []
accession_list = []
fo = open('species2.txt', 'w')
fi = open('species.txt', 'r')
fi.readline()
fo.write('Accession\tSpecies\tNr of reads\n')
for line in fi:
    accession = line.split('\t\t')[0]
    species = line.split('\t\t')[1]
    nor = line.split('\t\t')[2]
    if species in species_list:
        nor_list[species.index(species)] += int(nor)
    else:
        print(species)
        species_list.append(species)
        accession_list.append(accession)
        nor_list.append(int(nor))
for a, s, n in zip(accession_list, species_list, nor_list):
    fo.write(a + '\t' + s + '\t' + str(n) + '\n')
fi.close()
fo.close()
