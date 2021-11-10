import time
import datetime

start_time = time.time()
if __name__ == "__main__":
    import argparse
    from collections import OrderedDict
    PARSER = argparse.ArgumentParser(description='Assemble reads to species.')
    PARSER.add_argument('-f', type=str, help='read list file', required=True)
    PARSER.add_argument('-o', type=str, help='read list file', required=True)
    ARGS = PARSER.parse_args()
    with open(ARGS.f) as fi:
        inp_data = fi.read()
        reads_per_taxId = dict()
        totReads = 0
        for item in inp_data.split('\n'):
            try:
                nr_reads = int(item.split('\t')[0].split(';size=')[1])
                totReads += nr_reads
                taxId = item.split('\t')[1]
                if taxId in reads_per_taxId.keys():
                    reads_per_taxId[taxId] += nr_reads
                else:
                    reads_per_taxId[taxId] = nr_reads
            except:
                break
    ordered_taxIds = OrderedDict(sorted(reads_per_taxId.items(), key=lambda t: t[1], reverse=True))
    with open(ARGS.o, 'w') as fi:
        for key in ordered_taxIds:
            if ordered_taxIds[key] > 0:
                fi.write('{}: {} ({}%)\n'.format(key, ordered_taxIds[key], round(100*(ordered_taxIds[key]/totReads),
                                                                                 2)))
        fi.write('Total number of reads: {}'.format(totReads))
    time_in_seconds = time.time() - start_time
    # hours, minutes = divmod(time_in_seconds, 3600)
    # minutes, seconds = divmod(minutes, 60)
    print('Completed in {}'.format(str(datetime.timedelta(seconds=int(time_in_seconds)))))



