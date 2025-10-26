#!/usr/bin/env python3
import sys
import gzip
import re

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def rev_edge_info(edge):
    direction_dict = {"+" : "-", "-" : "+"}
    regex = r'(\d+)([-+])(\d+)([-+])'
    match = re.match(regex, edge)
    node1, direction1, node2, direction2 = match.groups()    
    rev_edge = str(node2) + direction_dict[str(direction2)] + str(node1) + direction_dict[str(direction1)]
    return rev_edge

def get_edge_seq(seq, direct):
    if direct == ">":
        return seq
    elif direct == "<":
        rev_dict = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
        rev_seq = ''.join([rev_dict[base] for base in seq[::-1]])
        return rev_seq

gfa_file = sys.argv[1]
target_sample = sys.argv[2]
target_hap = sys.argv[3]
out_fa_file = sys.argv[4]

node_seq_dict = {}
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        seq = line.strip().split()[2]
        node_seq_dict[node] = seq

fo = open(out_fa_file, "w")
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.strip().split()[1]
        hap = line.strip().split()[2]
        if sample == target_sample and hap == target_hap:
            contig = line.strip().split()[3]
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            seq = ""
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                seq += get_edge_seq(node_seq_dict[node], direct)
            fo.write(">" + contig + "\n")
            fo.write(seq + "\n")
fo.close()

