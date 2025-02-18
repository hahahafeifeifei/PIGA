#!/usr/bin/env python3
import sys
import gzip
import re
from collections import Counter


def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")
    

gfa_file = sys.argv[1]
max_sum_coverage = int(sys.argv[2])
max_hap_coverage = int(sys.argv[3])
out_bed_file = sys.argv[4]


node_len_dict = {}
node_coverage_dict = {}
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        node_coverage_dict[node] = 0
        node_len_dict[node] = len(line.strip().split()[2])


filter_node_dict = {}
pre_hap_name = ""
hap_node = []
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.strip().split()[1]
        if sample != "_MINIGRAPH_":
            hap = line.strip().split()[2]
            hap_name = sample + "#" + hap
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            if hap_name != pre_hap_name and pre_hap_name != "":
                hap_node_count = Counter(hap_node)
                for node, count in hap_node_count.items():
                    node_coverage_dict[node] += count
                    if count > max_hap_coverage:
                        filter_node_dict[node] = True
                hap_node = []
            pre_hap_name = hap_name
            hap_node += path_node
if sample != "_MINIGRAPH_":
    hap_node_count = Counter(hap_node)
    for node, count in hap_node_count.items():
        node_coverage_dict[node] += count
        if count > max_hap_coverage:
            filter_node_dict[node] = True

for node, coverage in node_coverage_dict.items():
    if coverage > max_sum_coverage:
        filter_node_dict[node] = True

fo = open(out_bed_file, "w")
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.strip().split()[1]
        contig = line.strip().split()[3]
        contig_start = int(line.strip().split()[4])
        path_info = line.strip().split()[6]
        path_node = re.split('<|>', path_info)[1:]

        pre_i = 0
        pre_pos = contig_start
        pos = contig_start
        for i in range(len(path_node)):
            node = path_node[i]
            if node not in filter_node_dict:
                if pre_i != i:
                    fo.write(contig + '\t' + str(pre_pos) + '\t' + str(pos) + '\n')
                pre_i = i + 1
                pre_pos = pos + node_len_dict[node]
            pos += node_len_dict[node]
        if pre_i != i + 1:
            fo.write(contig + '\t' + str(pre_pos) + '\t' + str(pos) + '\n')
fo.close()


