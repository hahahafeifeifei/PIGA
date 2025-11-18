#!/usr/bin/env python3

# Pacbio median reads extraction script
#
# Extract the median reads of each zmw from fasta file. The fasta file must be sorted by the zmw name.
# 
# Date: 2022/12/10

import sys
import gzip
import numpy as np

input_file = sys.argv[1]
output_file = sys.argv[2]

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def writefile(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, "wt")
    else:
        return open(filename, "w")

def choose_median_zmw(subreads_len_list):
    if len(subreads_len_list) <= 3:
        subreads_pos=np.where(subreads_len_list == np.amax(subreads_len_list))[0][0]
    else:
        median_distance = abs(subreads_len_list - np.median(subreads_len_list))
        median_pos = np.where(median_distance == np.amin(median_distance))[0]
        subreads_pos = np.where(subreads_len_list == np.amax(np.array(subreads_len_list)[median_pos]))[0][0]
    return subreads_pos

fo = writefile(output_file)

zmw_list = []
subreads_len_list = []
subreads_list = []
pre_zmw_name = "null"
line_number = 0

for line in openfile(input_file):
    line_number += 1
    if line_number%2 == 1:
        subreads_list.append(line.strip())
        subreads_name = line.strip().strip(">")
        zmw_name = subreads_name.split("/")[0]+"/"+subreads_name.split("/")[1]
        
        if pre_zmw_name != "null" and zmw_name != pre_zmw_name:
            pos = choose_median_zmw(subreads_len_list)
            zmw_reads = zmw_list[pos]  
            zmw_list = []
            subreads_len_list = []

            fo.write("\n".join(zmw_reads))
            fo.write("\n")
        
        pre_zmw_name = zmw_name

    elif line_number%2 == 0:
        subreads_len = len(line.strip())
        subreads_len_list.append(subreads_len)
        subreads_list.append(line.strip())
        zmw_list.append(subreads_list)
        subreads_list = []

pos = choose_median_zmw(subreads_len_list)
zmw_reads = zmw_list[pos]
fo.write("\n".join(zmw_reads))
fo.write("\n")
fo.close()
