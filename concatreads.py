import os
import argparse


def get_variants(sample, marker):
    variants = []
    variant = 'a'  # absent
    fraction = 1
    seq = ''
    for item in marker:
        if item[0].split('_')[0] == sample:
            variant = item[0].split('_')[1]
            if variant == 'f:100%':
                variant = 's'  # single
            # TODO Why split at '_c', no '_c' exist in input?
            fraction = float(item[0].split('f:')[1].split('_c')[0][:-1])/100
            seq = item[1]
            variants.append([sample, variant, fraction, seq])
    if not variants:
        variants = [[sample, variant, fraction, seq]]
    return variants


def get_sample_names():
    sample_names = set()
    for marker in markerlst:
        for sample in marker:
            sample_names.add(sample[0].split('_')[0])
    return sample_names


def get_concat_variants(sample):
    variants_old = [['', 1, '']]
    variants_new = []
    for marker in markerlst:
        variants = get_variants(sample, marker)
        for variant in variants:
            for item in variants_old:
                new_label = item[0] + variant[1]
                new_fraction = item[1] * variant[2]
                new_seq = item[2] + variant[3]
                variants_new.append([new_label, new_fraction, new_seq])
        variants_old = variants_new.copy()
        variants_new = []
    print(variants_old)
    test = input('wait')
    return variants_old


PARSER = argparse.ArgumentParser(description='concatenate markers')
PARSER.add_argument('-d', type=str, help='input file directory', required=True)
PARSER.add_argument('-f', type=str, help='output filename', required=True)
ARGS = PARSER.parse_args()
markerlst = list()
for filename in os.listdir(ARGS.d):
    if filename.endswith('.afa'):
        with open(ARGS.d + '\\' + filename) as fi:
            markerlst.append([x.split('\n')[:2] for x in fi.read().split('>')][1:])
outf = open(ARGS.f, 'w')
for sample in get_sample_names():
    if 'a' not in get_concat_variants(sample)[0][0] and not \
            get_concat_variants(sample)[0][0].count('v') > 1:
        for count, item in enumerate(get_concat_variants(sample)):
            result = '>' + sample
            if 'v' in item[0]:
                result += '_v' + str(count + 1) + '_' + str(round(item[1]*100)) + '%'
            result += '\n' + item[2] + '\n'
            outf.write(result)
outf.close()
