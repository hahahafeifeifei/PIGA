#!/usr/bin/env python3

# Shift pos of aligments and transform to bed file
#
# Extract the median reads of each zmw from fasta file. The fasta file must be sorted by the zmw name.
# 
# Date: 2022/12/11

import pysam
import sys

for read in pysam.AlignmentFile(sys.stdin, "r"):
    if (not read.is_secondary) and (not read.is_unmapped):
        ref_chr = read.reference_name.split(":")[0]
        ref_start = int(read.reference_name.split(":")[1].split("-")[0])
        read_start = ref_start + read.reference_start - 1
        read_end = ref_start + read.reference_end - 1
        read_name = read.qname
        print(ref_chr + '\t' + str(read_start) + '\t' + str(read_end) + '\t' + read_name)
