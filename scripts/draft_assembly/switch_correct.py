#!/usr/bin/env python3
import sys
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

vcf_file = sys.argv[1]
fo_vcf = open(sys.argv[2],"w")
fo_bed = open(sys.argv[3],"w")

pre_chr = "0"
pre_pos = 0

for line in openfile(vcf_file):
    if line[0:2] == "##":
        fo_vcf.write(line)
    elif line[0] == "#":
        fo_vcf.write('\t'.join(line.split()[:10]) + '\n')
    else:
        chr = line.split()[0]
        pos = line.split()[1]
        if pre_chr != chr:
            switch = 1
        
        gt1 = line.strip().split()[9].split(":")[0]
        gt2 = line.strip().split()[10].split(":")[0]
        if "|" not in gt1:
            continue
        if gt2 != "./.":
            if (switch == 1 and gt1 != gt2) or (switch == -1 and gt1 == gt2):
                switch *= -1
                fo_bed.write(chr + '\t' + str(0 if int(pre_pos)-25001 < 0 else int(pre_pos)-25001) + '\t' + str(int(pos)+25000) + '\n')
        
        if switch == 1:
            fo_vcf.write('\t'.join(line.split()[:8]) + '\t' + "GT" + "\t" + gt1 + '\n')
        elif switch == -1:
            rev_gt1 = gt1.split("|")[1] + '|' + gt1.split("|")[0]
            fo_vcf.write('\t'.join(line.split()[:8]) + '\t' + "GT" + "\t" + rev_gt1 + '\n')

        pre_chr = chr
        pre_pos = pos

fo_vcf.close()
fo_bed.close()
