#!/usr/bin/env python3
import sys
import gzip
import re
import subprocess

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
node_dict = {}
edge_dict = {}

#####load the node and edge data
for line in openfile(gfa_file):
    if line[0] == "W":
        path_info = line.strip().split()[6]
        path_node = re.split('<|>', path_info)[1:]
        path_direct = re.split('\d+', path_info)[:-1]
        pre_node_info = "null"
        for i in range(len(path_node)):
            node = path_node[i]
            direct = path_direct[i]
            node_info = node + direct.replace(">","+").replace("<","-")
            node_dict[node] = True
            if pre_node_info != "null":
                pre_node = pre_node_info[:-1]
                edge = pre_node_info + node_info
                edge_dict[edge] = True
            pre_node_info = node_info

fo = open(sys.argv[2], "w")
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        if node in node_dict:
            fo.write(line)
            del node_dict[node]
    elif line[0] == "L":
        edge_info = line.strip().split()
        edge = edge_info[1] + edge_info[2] + edge_info[3] + edge_info[4]
        rev_edge = rev_edge_info(edge)
        if edge in edge_dict:
            fo.write(line)
            del edge_dict[edge]
        elif rev_edge in edge_dict:
            fo.write(line)
            del edge_dict[rev_edge]
    elif line[0] == "W":
        fo.write(line)
fo.close()
