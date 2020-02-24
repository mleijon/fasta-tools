import os


def get_variants(sample, marker):
    variants = []
    variant = 'a'
    fraction = 1
    seq = ''
    for item in marker:
        if item[0].split('_')[0] == sample:
            variant = item[0].split('_')[1]
            if variant == 'f:100%':
                variant = 's'
            fraction = float(item[0].split('f:')[1][:-1])/100
            seq = item[1]
            variants.append([sample, variant, fraction, seq])
    if not variants:
        variants = [[sample, variant, fraction, seq]]

    return variants


markerlst = []
variants_old = [['', 1, '']]
variants_new = []
for filename in os.listdir("C:\cp"):
    with open(filename) as fi:
        markerlst.append([x.split('\n')[:2] for x in fi.read().split('>')][1:])

for marker in markerlst:
    variants = get_variants('Cp.01.02.04', marker)
    for variant in variants:
        for item in variants_old:
            new_label = item[0] + variant[1]
            new_fraction = item[1] * variant[2]
            new_seq = item[2] + variant[3]
            variants_new.append([new_label, new_fraction, new_seq])
    variants_old = variants_new.copy()
    variants_new = []


