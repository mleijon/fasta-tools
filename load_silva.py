import gzip, urllib.request, zipfile, io, shutil, os

surl = "https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_" \
       "SSURef_NR99_tax_silva.fasta.gz"
turl = "https://ftp.arb-silva.de/release_138.1/Exports/taxonomy/taxmap" \
       "_slv_ssu_ref_nr_138.1.txt.gz"
nurl = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
print("Downloading " + nurl[nurl.rfind("/") + 1:] + " may take some time ... ",
      end="", flush=True)
allowedRanks = {"superkingdom": "k__", "phylum": "p__", "class": "c__",
                "order": "o__", "family": "f__", "genus": "g__",
                "species": "s__"}


def sp(line):
    return line.replace(b"\n", b"\t").split(b"\t|\t")


with zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(nurl).read()))\
        as zip_ref:
    with zip_ref.open([name for name in zip_ref.namelist() if
                       os.path.basename(name) == "nodes.dmp"][0]) as zf:
        nodes = {sp(line)[0]: [sp(line)[1], sp(line)[2].decode("UTF-8"), ""]
                 for line in zf}
    with zip_ref.open([name for name in zip_ref.namelist() if
                       name.endswith("names.dmp")][0]) as zf:
        for line in zf:
            s = sp(line)
            if s[3] == b"scientific name":
                nodes[s[0]][2] = s[1].decode("UTF-8")


def getLineage(byteTaxID):
    lin = {r: v for r, v in allowedRanks.items()}
    pid = byteTaxID
    if pid in nodes:
        mid = nodes[pid]
        while pid != b"1" and pid != mid[0]:
            if mid[1] in allowedRanks:
                lin[mid[1]] += mid[2]
            pid = mid[0]
            mid = nodes[pid]
    return "; ".join(v for k, v in lin.items())


print("done")
oname1 = surl[surl.rfind("/") + 1:].replace("fasta.gz", "fa.gz")
oname2 = oname1.replace("fa.gz", "txt")
print("Downloading " + turl[turl.rfind("/") + 1:] +
      " may take some time ... ", end="", flush=True)
with gzip.GzipFile(fileobj=urllib.request.urlopen(turl)) as gzTax,\
        open(oname2, 'w') as tO:
    next(gzTax)
    for line in gzTax:
        sp = line.strip().split(b"\t")
        tO.write(sp[0].decode("UTF-8") + "." + sp[1].decode("UTF-8") + "." +
                 sp[2].decode("UTF-8") + "\t" + getLineage(
            sp[5]) + "\n")
print("done")
print("Taxonomy output: " + oname2)
print("Downloading " + surl[surl.rfind("/") + 1:] +
      " may take some time ... ", end="", flush=True)
with gzip.GzipFile(fileobj=urllib.request.urlopen(surl)) as\
        gzSilva, gzip.open(oname1, 'wb') as fO:
    for line in gzSilva:
        if line.startswith(b">"):
            fO.write(line[:line.rfind(b" ", 0, line.find(b";"))] + b"\n")
        else:
            fO.write(line.replace(b"U", b"T"))
print("done")
print("Fasta output:    " + oname1)
