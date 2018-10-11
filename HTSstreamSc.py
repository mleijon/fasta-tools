#!/usr/bin/python3

from glob import glob
import os
cleaning = open("cleaning_commands.sh", 'w')
for r1 in glob("./00-rawdata/*_R1_*.gz"):
    r2 = r1.replace("R1", "R2")
    s = r1.split('/')[-1].replace("_L001_R1_001.fastq.gz", '')
    cmd = "hts_SuperDeduper -L ./01-Cleaned/" + s + "_stats.log -1 " + r1 + " -2 " + r2 + " -O | "
    cmd += "hts_SeqScreener -S -O -A -L ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_PolyATTrim -m 100 -S -O -A -L ./01-Cleaned/" + s + "_stats.log | "
    cmd += "hts_Stats -N phix-remover-adapters -S -A -L ./01-Cleaned/" + s + "_StatsStats.log -g -p ./01-Cleaned/" + s
    cleaning.write(cmd+'\n')
cleaning.close()

#      cmd += "hts_SeqScreener -S -O -A -L ./01-Cleaned/" + s + "_stats.log --seq adapters.fa -k 15 -x .01 | "
