#!/usr/bin/python

if __name__ == "__main__":
    import sys
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)
    fasta_in = sys.stdin.read()
    indata = fasta_in.replace(';', ':0;').split('>')
    del indata[0]
    indata_alt = []
    for i in range(len(indata)):
        indata_alt.append(indata[i].replace('NNNNNNNN', '\n@' + indata[i].split('\n')[0].replace(':0;', ':1;') + '@'))
    indata = []
    for item in indata_alt:
        part1 = '>' + item.split('\n', 1)[0]
        part2 = item.split('\n', 1)[1].replace('\n', '').split('@')
        new_item = part1 + '\n' + part2[0] + '\n' + '>' + part2[1] + '\n' + part2[2] + '\n'
        sys.stdout.write(new_item)
