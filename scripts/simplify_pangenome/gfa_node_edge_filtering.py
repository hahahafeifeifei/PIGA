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

gfa_file = sys.argv[1]
node_file = sys.argv[2]
edge_file = sys.argv[3]
ref_prefix = sys.argv[4]
out_gfa_file = sys.argv[5]

filter_node_dict = {}
filter_edge_dict = {}
for line in openfile(node_file):
    filter_node_dict[line.strip()] = True
for line in openfile(edge_file):
    filter_edge_dict[line.strip()] = True

node_len_dict = {}
ref_list = ref_prefix.split(",")
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        node_len_dict[node] = len(line.strip().split()[2])
    if line[0] == "W":
        sample = line.split()[1]
        if sample in ref_list:
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_node_info = "null"
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                node_info = node + direct.replace(">","+").replace("<","-")
                if node in filter_node_dict:
                    del filter_node_dict[node]
                if pre_node_info != "null":
                    edge = pre_node_info + node_info
                    rev_edge = rev_edge_info(edge)
                    if edge in filter_edge_dict:
                        del filter_edge_dict[edge]
                    elif rev_edge in filter_edge_dict:
                        del filter_edge_dict[rev_edge]
                pre_node_info = node_info

fo = open(out_gfa_file, "w")
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        if node not in filter_node_dict:
            fo.write(line)
    if line[0] == "L":
        edge_info = line.strip().split()
        edge = edge_info[1] + edge_info[2] + edge_info[3] + edge_info[4]
        if edge_info[1] not in filter_node_dict and edge_info[3] not in filter_node_dict and edge not in filter_edge_dict:
            fo.write(line)
    if line[0] == "W":
        sample = line.strip().split()[1]
        hap = line.strip().split()[2]
        contig = line.strip().split()[3]
        contig_start = int(line.strip().split()[4])
        contig_end = int(line.strip().split()[5])
        path_info = line.strip().split()[6]

        path_node = re.split('<|>', path_info)[1:]
        path_direct = re.split('\d+', path_info)[:-1]
        pre_node_info = "null"
        pre_i = 0
        pre_pos = contig_start
        pos = contig_start
        for i in range(len(path_node)):
            node = path_node[i]
            direct = path_direct[i]            
            node_info = node + direct.replace(">","+").replace("<","-")

            if pre_node_info != "null":
                edge = pre_node_info + node_info
                rev_edge = rev_edge_info(edge)
                if edge in filter_edge_dict or rev_edge in filter_edge_dict:
                    if pre_i != i:
                        write_path = ''.join([path_direct[j] + path_node[j] for j in range(pre_i, i)])
                        fo.write('\t'.join(["W", sample, hap, contig, str(pre_pos), str(pos), write_path]) + '\n')
                    pre_i = i
                    pre_pos = pos
            if node in filter_node_dict:
                if pre_i != i:
                    write_path = ''.join([path_direct[j] + path_node[j] for j in range(pre_i, i)])
                    fo.write('\t'.join(["W", sample, hap, contig, str(pre_pos), str(pos), write_path]) + '\n')
                pre_i = i + 1
                pre_pos = pos + node_len_dict[node]
            pos += node_len_dict[node]
            pre_node_info = node_info

        if pre_i != i + 1:
            write_path = ''.join([path_direct[j] + path_node[j]for j in range(pre_i, i + 1)])
            fo.write('\t'.join(["W", sample, hap, contig, str(pre_pos), str(pos), write_path]) + '\n')

fo.close()

