#!/usr/bin/python3
"""Demultiplex NGS data with multiplexed PCR amplicons listed in a fasta-file
"""
from fasta import FastaList


def sep(opsys):
    if opsys == 'nt':
        return '\\'
    else:
        return '/'


if __name__ == "__main__":
    import copy
    import collections
    import argparse
    import shutil
    import os
    import sys
    import operator

    PARSER = argparse.ArgumentParser(description='Removes redundat sequences'
                                                 'after manual trimming of '
                                                 'the result of demultiplex.py'
                                                 'and updates count.')
    PARSER.add_argument('-id', type=str, help='Directory for input data',
                        required=True)
    PARSER.add_argument('-od', type=str, help='Directory for output data',
                        required=True)
    PARSER.add_argument('-c', action='store_true',
                        help='if set counts will be shown', default=False)
    PARSER.add_argument('-r', action='store_true',
                        help='if set reference will be shown', default=False)
    PARSER.add_argument('-t', type=float, help='threshold for fraction of the'
                                               'total amounts of reads for a'
                                               'sample', default=0)

    ARGS = PARSER.parse_args()
    # Some controls of input data
    if not ((ARGS.id.endswith(sep(os.name)) and
             ARGS.od.endswith(sep(os.name)))):
        sys.exit('Invalid directory name. Exits.')
    if os.path.isfile(ARGS.od[:-1]):
        print('{} is a file'.format(ARGS.od[:-1]))
        sys.exit('Exits')
    if os.path.isdir(ARGS.od):
        shutil.rmtree(ARGS.od)
    os.mkdir(ARGS.od)
    # loop over files produced by demultiplex.py and with manually removed indels
    # and cropped to identical length
    for seqfile in os.listdir(ARGS.id):
        if seqfile.endswith('.afa'):

            # Dict with sequence_names (sequences names are sample names with a
            # variant incremental number (#i) attached (Sample_name_#i)) as keys and
            # a list [sample_name, sequence, sequnce count] as values. These are
            # read from the fasta files produced by demultiplex.py
            sample = dict()

            # Dict with sample_names as keys containg  {sample (dict()) items} (i.e.
            # sequence variants)
            sample_reduced = dict(dict())

            inp_seqs = FastaList(ARGS.id + seqfile)
            for seq in inp_seqs.seq_list:

                # Name with incremental number (_#i)
                seq_name = seq.split('_count:')[0][1:]
                if 'ref' in seq_name:
                    seq_name = seq.split('\n')[0][1:]
                    seq_count = 1
                else:
                    seq_count = int(seq.split('_count:')[1].split('_')[0])
                dna_seq = seq.split('\n')[1]

                # Name without incremental number
                if 'ref' in seq_name:
                    sample_name = seq_name
                else:
                    sample_name = seq_name.rsplit('_', 1)[0]
                sample[seq_name] = [sample_name, dna_seq, seq_count]

            # Loop over all sequences of the input fasta file
            for seqname in sample:
                current_sample = sample[seqname][0]
                found = False
                if current_sample in sample_reduced.keys():

                    # Loop over all sequences collected for a sample (samplename)
                    # in Dict sample_reduced if the sequence exist add the counts
                    # if the sequence is new add it to the dictionary
                    for seq_variant in sample_reduced[current_sample]:
                        if sample[seqname][1] ==\
                                sample_reduced[current_sample][seq_variant][1]:
                            sample_reduced[current_sample][seq_variant][2] +=\
                                sample[seqname][2]
                            found = True
                    if not found:
                        sample_reduced[current_sample][seqname] = sample[seqname]
                else:
                    sample_reduced[current_sample] = {seqname: sample[seqname]}

            # Rename the sample_names (keys of sample_reduced) to make sure sorting
            # will work. Sort make a sorted version of sample_reduced:
            # samples_red_ord
            sample_tmp = dict()
            for key in sample_reduced:
                if 'ref' in key:
                    name = key
                else:
                    n = key.split('_')
                    name = n[0] + '_' + n[1].zfill(2) + '_' + n[2].\
                        zfill(2) + '_' + n[3].zfill(2)
                sample_tmp[name] = sample_reduced[key]
            sample_reduced = copy.deepcopy(sample_tmp)
            samples_red_ord = collections.OrderedDict(sorted(
                sample_reduced.items()))

            # Carries the variants of sample_name in a nested list with the members:
            # [variant_name, variant_sequence, variant_count]. The list is sorted in
            # reverse (highest first) by the variant_count.
            result = []
            tmp_result = []

            # Ordered (by sample_name) triple netsted list with all results:
            # [[[variant1], [variant2],...], [[variant1], [variant2]]]. The variants
            # for each sample are ordered by count
            result_all = []

            # Create a new nested list 'new_result' from 'result' that adds the
            # total count_sum for the variants to the variants listed in new_result:
            # [variant_name, variant_sequence, variant_count, count_sum] and
            # modifies the variants name by removing the incremental number and add
            # a _v#, where # is a version number for the variant. the nested lists
            # of 'new_result' are appended to the final result list: 'result_all'
            for item in samples_red_ord:
                count_sum = 0
                for var in samples_red_ord[item]:
                    count_sum += samples_red_ord[item][var][2]
                    tmp_result.append([var, samples_red_ord[item][var][1],
                                       samples_red_ord[item][var][2]])
                    tmp_result.sort(key=operator.itemgetter(2), reverse=True)

                # Exclude variant below the threshold set by -t (default =0). Adjust
                # the count_sum accordingly.
                for var in tmp_result:
                    if (var[2]/count_sum) < ARGS.t:
                        count_sum -= var[2]
                    else:
                        result.append(var)

                new_result = []
                variant_count = 0
                for variant in result:
                    variant.append(count_sum)
                    if len(result) > 1:
                        variant_count += 1
                        variant[0] = variant[0].rsplit('_', 1)[0] + '_v' + str(
                            variant_count)
                    else:
                        if 'ref' not in variant[0]:
                            variant[0] = variant[0].rsplit('_', 1)[0]
                        else:
                            pass
                    new_result.append(variant)
                result_all.append(new_result)
                result = []
                tmp_result = []

            with open(ARGS.od + seqfile, 'w') as fi:
                for item1 in result_all:
                    for item2 in item1:
                        if 'ref' in item2[0]:
                            name = item2[0]
                            n = [name]
                        else:
                            n = item2[0].split('_')
                            name = n[0] + '.' + n[1].zfill(2) + '.' + n[2].\
                                zfill(2) + '.' + n[3].zfill(2)
                        if len(n) > 4:  # variants exit
                            name += '_' + n[4]
                        dna_seq = item2[1]
                        seq_count = item2[2]
                        seq_count_sum = item2[3]
                        if 'ref' in name and not ARGS.r:
                            continue
                        if ARGS.c:
                            fi.write('>{}_f:{}%_c:{}\n{}\n'.format(name, round(100*(
                                    seq_count/seq_count_sum)), seq_count, dna_seq))
                        else:
                            fi.write('>{}_f:{}%\n{}\n'.format(name, round(100*(
                                    seq_count/seq_count_sum)), dna_seq))
