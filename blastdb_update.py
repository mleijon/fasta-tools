import argparse
import subprocess
from os.path import exists as file_exists

url_tax_aa = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz'
url_tax_nt = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
url_gb_release = 'https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number'
args = ['wget', '--quiet', '--timestamping']


def md5_sum(file_cont):
    import hashlib
    md5_hash = hashlib.md5()
    md5_hash.update(file_cont)
    result = md5_hash.hexdigest()
    return result


def get_file(url, md5):
    global args
    args.append(url)
    file_name = args[3].rsplit('/', maxsplit=1)[1]
    print('Working with file: {}'.format(file_name))
    subprocess.call(args)
    if md5:
        with open(file_name, 'rb') as f:
            f_md5sum = md5_sum(f.read())
        args.append(args.pop() + '.md5')
        subprocess.call(args)
        with open(file_name + '.md5') as f:
            s_md5sum = f.read().split()[0]
        if s_md5sum == f_md5sum:
            print("Checksum OK!")
        else:
            print('Checksum error!')
    else:
        print('No checksum file.')
    args.pop()


def makeblastdb(file_str):
    args = ['makeblastdb -in ' + file_str + ' -input_type asn1_bin -dbtype nucl -parse_seqids -out test']
    subprocess.call(args, shell=True)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-p', type=str, help='input genbank patition', required=True)
    PARSER.add_argument('-t', type=str, help='type of db (nt/aa)', required=True)
    ARGS = PARSER.parse_args()
    filestr = ''
    url_gen = str
    get_file(url_gb_release, False)
    print('** PHASE 1 - Checking sequence files **')
    for nr in range(1, 3):
        if ARGS.t.lower() == 'prot':
            try:
                url_gen = 'https://ftp.ncbi.nlm.nih.gov/ncbi-asn1/protein_fasta/gb' + ARGS.p + str(nr) + '.fsa_aa.gz'
                args.append(url_gen)
                subprocess.check_call(args)
            except subprocess.CalledProcessError:
                break
        elif ARGS.t.lower() == 'nucl':
            try:
                url_gen = 'https://ftp.ncbi.nlm.nih.gov/ncbi-asn1/gb' + ARGS.p + str(nr) + '.aso.gz'
                args.append(url_gen)
                subprocess.check_call(args)
            except subprocess.CalledProcessError:
                break
        print(args[3].rsplit('/', maxsplit=1)[1])
        filestr += url_gen.rsplit('/', maxsplit=1)[1][:-3] + ' '
        args.pop()
    filestr = '"' + filestr[:-1] + '"'
    print('\nDone checking sequence files!\nFiles updated if newer versions were found on server.\n')
    print('** PHASE 2 - Checking taxonomy file **')
    get_file(url_tax_aa, True)
    get_file(url_tax_nt, True)
    makeblastdb(filestr)


