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

def get_edge_seq(seq, direct):
    if direct == "+":
        return seq
    elif direct == "-":
        rev_dict = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
        rev_seq = ''.join([rev_dict[base] for base in seq[::-1]])
        return rev_seq

gfa_file = sys.argv[1]
ref_samples = sys.argv[2].split(",")
minigraph_sample = sys.argv[3].split(",")
tmp_dir = sys.argv[4]

node_seq_dict = {}
node_homopolymer_dict = {}
node_dict = {}
edge_dict = {}
node_sample_dict = {}
edge_sample_dict = {}

#####load the node and edge data
print("1. Load the graph information ...")
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        seq = line.strip().split()[2].upper()
        node_dict[node] = [[]] * 8
        node_dict[node][1] = 0
        node_dict[node][2] = 0
        node_dict[node][3] = 0
        node_dict[node][4] = 0
        node_dict[node][5] = 0
        node_dict[node][6] = len(line.strip().split()[2])
        node_homopolymer_dict[node] = [0, 0]
        node_seq_dict[node] = seq
    elif line[0] == "L":
        edge_info = line.strip().split()
        edge = edge_info[1] + edge_info[2] + edge_info[3] + edge_info[4]
        edge_dict[edge] = [[]] * 8
        edge_dict[edge][2] = 0
        edge_dict[edge][3] = 0
        edge_dict[edge][4] = 0
        edge_dict[edge][5] = 0
        edge_dict[edge][6] = 0

#####detect the homopolymer edge
for line in openfile(gfa_file):
    if line[0] == "L":
        edge_info = line.strip().split()
        edge = edge_info[1] + edge_info[2] + edge_info[3] + edge_info[4]
        node1 = edge_info[1]
        node2 = edge_info[3]
        direct1 = edge_info[2]
        direct2 = edge_info[4]
        node1_seq = get_edge_seq(node_seq_dict[node1], direct1)
        node2_seq = get_edge_seq(node_seq_dict[node2], direct2)
        if node1_seq[-1] == node2_seq[0] and node2_seq[0] != "N":
            edge_dict[edge][7] = 1
        else:
            edge_dict[edge][7] = 0


####extract the graphlet feature
print("2. Run ORCA to extract the graphlet (motif) features ...")
i = 0
node_index_dict = {}
for node in node_dict:
    node_index_dict[node] = i
    i += 1

edge_orca_dict = {}
for edge in edge_dict:
    edge_node = re.split(r'[+-]', edge)[:-1]
    node1 = edge_node[0]
    node2 = edge_node[1]
    if int(node1) >= int(node2):
        tmp = node2
        node2 = node1
        node1 = tmp
    index1 = node_index_dict[node1]
    index2 = node_index_dict[node2]
    if index1 != index2:
        edge_index = str(index1) + ' ' + str(index2)
        if edge_index not in edge_orca_dict:
            edge_orca_dict[edge_index] = [edge]
        else:
            edge_orca_dict[edge_index].append(edge)

fo = open(tmp_dir + "/tmp.orca.edge","w")
fo.write(str(len(node_index_dict)) + ' ' + str(len(edge_orca_dict)) + '\n')
for edge_index in edge_orca_dict:
    index1 = edge_index.split()[0]
    index2 = edge_index.split()[1]
    fo.write(index1 + ' ' + index2 + '\n')
fo.close()

subprocess.run("orca.exe node 5 {} {}".format(tmp_dir + "/tmp.orca.edge", tmp_dir + "/tmp.orca.node.graphlet"), shell = True, stdout = None, stderr = None)
subprocess.run("orca.exe edge 5 {} {}".format(tmp_dir + "/tmp.orca.edge", tmp_dir + "/tmp.orca.edge.graphlet"), shell = True, stdout = None, stderr = None)
subprocess.run("rm {}".format(tmp_dir + "/tmp.orca.edge"), shell = True, stdout = None, stderr = None)

node_list = list(node_index_dict.keys())
i = 0
for line in openfile(tmp_dir + "/tmp.orca.node.graphlet"):
    node = node_list[i]
    node_dict[node][0] = node_index_dict[node]
    graphlet_list = line.strip().split()
    for graphlet in graphlet_list:
        node_dict[node].append(graphlet)
    if len(graphlet_list) < 73:
        node_dict[node] += ['nan']*(73-len(graphlet_list))
    i += 1

edge_info_list = [[edge_list, edge_index.split()] for edge_index, edge_list in edge_orca_dict.items()]
i = 0
for line in openfile(tmp_dir + "/tmp.orca.edge.graphlet"):
    edge_list = edge_info_list[i][0]
    for edge in edge_list:
        edge_dict[edge][0] = edge_info_list[i][1][0]
        edge_dict[edge][1] = edge_info_list[i][1][1]
        graphlet_list = line.strip().split()
        for graphlet in graphlet_list:
            edge_dict[edge].append(graphlet)
        if len(graphlet_list) < 68:
            edge_dict[edge] += ['nan']*(68-len(graphlet_list))
    i += 1
subprocess.run("rm {} {}".format(tmp_dir + "/tmp.orca.node.graphlet", tmp_dir + "/tmp.orca.edge.graphlet"), shell = True, stdout = None, stderr = None)

#####extract the feature related to allele path
print("3. Extract the allele count related features ...")

pre_sample = "null"
node_sample_dict = {}
edge_sample_dict = {}
ref_node_dict = {}
ref_edge_dict = {}
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.split()[1]
        if sample != pre_sample: 
            if pre_sample not in ref_samples and pre_sample not in minigraph_sample and pre_sample != "null":
                for node in node_sample_dict:
                    ac = len(node_sample_dict[node])
                    if ac == 1:
                        node_dict[node][1] += 1
                        node_dict[node][3] += 1
                        node_dict[node][4] += 1
                    elif ac == 2:
                        node_dict[node][2] += 1
                        node_dict[node][3] += 2
                        node_dict[node][4] += 1
                node_sample_dict = {}
                for edge in edge_sample_dict:
                    ac = len(edge_sample_dict[edge])
                    if ac == 1:
                        edge_dict[edge][2] += 1
                        edge_dict[edge][4] += 1
                        edge_dict[edge][5] += 1
                    elif ac == 2:
                        edge_dict[edge][3] += 1
                        edge_dict[edge][4] += 2
                        edge_dict[edge][5] += 1
                edge_sample_dict = {}

        if sample in ref_samples:
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_node_info = "null"
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                node_info = node + direct.replace(">","+").replace("<","-")
                ref_node_dict[node] = True
                if pre_node_info != "null":
                    edge = pre_node_info + node_info
                    rev_edge = rev_edge_info(edge)
                    if edge in edge_dict:
                        ref_edge_dict[edge] = True
                    elif rev_edge in edge_dict:
                        ref_edge_dict[rev_edge] = True
                pre_node_info = node_info
        elif sample in minigraph_sample:
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_node_info = "null"
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                node_info = node + direct.replace(">","+").replace("<","-")
                node_dict[node][5] = 1
                if pre_node_info != "null":
                    edge = pre_node_info + node_info
                    rev_edge = rev_edge_info(edge)
                    if edge in edge_dict:
                        edge_dict[edge][6] = 1
                    elif rev_edge in edge_dict:
                        edge_dict[rev_edge][6] = 1
                pre_node_info = node_info
        else:
            hap = line.split()[2]
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_node_info = "null"
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                node_info = node + direct.replace(">","+").replace("<","-")
                if node not in node_sample_dict:
                    node_sample_dict[node] = {}
                node_sample_dict[node][hap] = True
                if pre_node_info != "null":
                    pre_node = pre_node_info[:-1]
                    edge = pre_node_info + node_info
                    rev_edge = rev_edge_info(edge)
                    if edge in edge_dict:
                        edge = edge
                    elif rev_edge in edge_dict:
                        edge = rev_edge  
                    if edge_dict[edge][7] == 1:
                        node_homopolymer_dict[node][0] += 1
                        node_homopolymer_dict[node][1] += 1
                        node_homopolymer_dict[pre_node][0] += 1
                        node_homopolymer_dict[pre_node][1] += 1
                    else:
                        node_homopolymer_dict[node][1] += 1
                        node_homopolymer_dict[pre_node][1] += 1
                    if edge not in edge_sample_dict:
                        edge_sample_dict[edge] = {}
                    edge_sample_dict[edge][hap] = True
                pre_node_info = node_info
        pre_sample = sample
    
if pre_sample not in ref_samples and pre_sample not in minigraph_sample and pre_sample != "null":
    for node in node_sample_dict:
        ac = len(node_sample_dict[node])
        if ac == 1:
            node_dict[node][1] += 1
            node_dict[node][3] += 1
            node_dict[node][4] += 1
        elif ac == 2:
            node_dict[node][2] += 1
            node_dict[node][3] += 2
            node_dict[node][4] += 1
    node_sample_dict = {}
    for edge in edge_sample_dict:
        ac = len(edge_sample_dict[edge])
        if ac == 1:
            edge_dict[edge][2] += 1
            edge_dict[edge][4] += 1
            edge_dict[edge][5] += 1
        elif ac == 2:
            edge_dict[edge][3] += 1
            edge_dict[edge][4] += 2
            edge_dict[edge][5] += 1
    edge_sample_dict = {}

print("4. Summarize and print features ...")
node_feature_file = open(sys.argv[5],"w")
for node in node_dict:
    if node in ref_node_dict:
        continue
    homopolymer = 0
    if node_homopolymer_dict[node][1] != 0:
        homopolymer = node_homopolymer_dict[node][0]/node_homopolymer_dict[node][1]
    homopolymer = round(homopolymer, 4)
    node_dict[node][7] = homopolymer
    node_feature_file.write("\t".join([node] + [str(info) for info in node_dict[node]]) + '\n')
node_feature_file.close()

edge_feature_file = open(sys.argv[6],"w")
for edge in edge_dict:
    if edge in ref_edge_dict:
        continue
    edge_node = re.split(r'[+-]', edge)[:-1]
    node1 = edge_node[0]
    node2 = edge_node[1]
    if node1 == node2:
        edge_dict[edge] += ['nan'] * 68
    edge_feature_file.write("\t".join([edge] + [str(info) for info in edge_dict[edge]]) + '\n')
edge_feature_file.close()
