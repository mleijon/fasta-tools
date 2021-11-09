import argparse
from fasta import FastaList
if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input fasta file',
                        required=True)
    PARSER.add_argument('--start', type=int, help='start cut', required=True)
    PARSER.add_argument('--end', type=int, help='end cut', required=True)
    ARGS = PARSER.parse_args()
    cut = ':' + str(ARGS.start) + '-' + str(ARGS.end)
    title = FastaList(ARGS.f).seq_list[0].split()[0] + cut + '\n'
    seq = FastaList(ARGS.f).seq_list[0].split()[1][ARGS.start:ARGS.end]
    outfile = ARGS.f.rsplit('.')[0] + cut.replace(':', '_') + '.' + ARGS.f.rsplit('.')[1]
    with open(outfile, 'w') as f:
        f.write(title)
        f.write(seq + '\n')
