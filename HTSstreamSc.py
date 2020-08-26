#!/usr/bin/python3

from glob import glob
cleaning = open("cleaning_commands.sh", 'w')
for r1 in glob("./00-rawdata/*_R1_*.gz"):
    r2 = r1.replace("R1", "R2")
    s = r1.split('/')[-1].replace("_L001_R1_001.fastq.gz", '')
    cmd = "hts_SuperDeduper -l 30 -L ./01-Cleaned/" + s + "_stats.log -1 " + r1 + " -2 " + r2 + " | "
    cmd += "hts_SeqScreener -A ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_PolyATTrim -M 100 -A ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_QWindowTrim -q 13 -A ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_Stats -N phix-remover-adapters -A ./01-Cleaned/" + s + "_stats.log -f ./01-Cleaned/" + s
    cleaning.write(cmd+'\n')
cleaning.close()