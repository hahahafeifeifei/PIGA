#!/usr/bin/env python3
import sys
import gzip
import re

#Identify the regions that are not aligned by minigraph path
#python3 gfa_nominigraph_bed.py input.gfa output.bed

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

gfa_file = sys.argv[1]
out_region_file = sys.argv[2]
filter_len = 100000

node_count_dict = {}
node_len_dict = {}
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        node_count_dict[node] = 0
        node_len_dict[node] = len(line.strip().split()[2])

mini_node_dict = {}
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.strip().split()[1]
        if sample == "_MINIGRAPH_":
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            for node in path_node:
                mini_node_dict[node] = True

fo = open(out_region_file, "w")
for line in openfile(gfa_file):
    if line[0] == "W":
        contig = line.strip().split()[3]
        path_info = line.strip().split()[6]
        path_node = re.split('<|>', path_info)[1:]
        pos = int(line.strip().split()[4])
        is_filter = False
        for node in path_node:
            if not is_filter:
                if node not in mini_node_dict:
                    is_filter = True
                    start_pos = pos
            elif is_filter:
                if node in mini_node_dict:
                    is_filter = False
                    end_pos = pos
                    if end_pos - start_pos > filter_len:
                        fo.write(contig + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n" )
            pos += node_len_dict[node]
        if is_filter:
            end_pos = pos
            if end_pos - start_pos > filter_len:
                fo.write(contig + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n" )
fo.close()