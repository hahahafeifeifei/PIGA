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
query_sample = sys.argv[2]
truth_sample = sys.argv[3]

node_dict = {}
edge_dict = {}
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        node_dict[node] = True
    elif line[0] == "L":
        edge_info = line.strip().split()
        edge = edge_info[1] + edge_info[2] + edge_info[3] + edge_info[4]
        edge_dict[edge] = True


query_node_dict = {}
query_edge_dict = {}
truth_node_dict = {}
truth_edge_dict = {}

for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.split()[1]
        if sample == query_sample:
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_node_info = "null"
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                node_info = node + direct.replace(">","+").replace("<","-")
                query_node_dict[node] = True
                if pre_node_info != "null":
                    pre_node = pre_node_info[:-1]
                    edge = pre_node_info + node_info
                    rev_edge = rev_edge_info(edge)
                    if edge in edge_dict:
                        query_edge_dict[edge] = True
                    elif rev_edge in edge_dict:
                        query_edge_dict[rev_edge] = True
                pre_node_info = node_info

        if sample == truth_sample:
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_node_info = "null"
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                node_info = node + direct.replace(">","+").replace("<","-")
                truth_node_dict[node] = True
                if pre_node_info != "null":
                    pre_node = pre_node_info[:-1]
                    edge = pre_node_info + node_info
                    rev_edge = rev_edge_info(edge)
                    if edge in edge_dict:
                        truth_edge_dict[edge] = True
                    elif rev_edge in edge_dict:
                        truth_edge_dict[rev_edge] = True
                pre_node_info = node_info


truth_node_file = open(sys.argv[4], "w")
for node in query_node_dict:
    if node in truth_node_dict:
        truth_node_file.write(node + "\t" + "TP" + "\n")
    else:
        truth_node_file.write(node + "\t" + "FP" + "\n")
truth_node_file.close()

truth_edge_file = open(sys.argv[5], "w")
for edge in query_edge_dict:
    if edge in truth_edge_dict:
        truth_edge_file.write(edge + "\t" + "TP" + "\n")
    else:
        truth_edge_file.write(edge + "\t" + "FP" + "\n")  
truth_edge_file.close()

