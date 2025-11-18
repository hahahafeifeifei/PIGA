#!/usr/bin/env python3

# Homopolymer extract script
#
# Extract the variants located in homopolymer region from vcf
# 
# Date: 2022/12/10

import sys

def homopolymer_filter(seq):
    compressed_seq = seq[0]
    for i in range(1,len(seq)):
        if seq[i] != seq[i-1]:
            compressed_seq += seq[i]
    return compressed_seq

for line in sys.stdin:
    if line[0] == "#":
        print(line.strip())
    else:
        ref = line.split('\t')[3].upper()
        alt = line.split('\t')[4].upper()
        compress_ref = homopolymer_filter(ref)
        compress_alt = homopolymer_filter(alt)
        if (compress_ref == compress_alt):
            print(line.strip())
