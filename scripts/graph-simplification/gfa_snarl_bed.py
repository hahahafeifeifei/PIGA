#!/usr/bin/env python3
import sys
import gzip
import re
import json

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")
    

gfa_file = sys.argv[1]
snarl_stat_file = sys.argv[2]
snarl_json_file = sys.argv[3]
min_size = int(sys.argv[4])
boundary_size = int(sys.argv[5])
out_snarl_bed_file = sys.argv[6]

snarl_list = []
snarl_son_dict = {}
for line in openfile(snarl_stat_file):
    line_info = line.strip().split()
    if line_info[0] != "Start":
        if int(line_info[6]) > min_size and int(line_info[7]) > min_size:
            snarl_id = '_'.join([str(x) for x in sorted([int(line_info[0]), int(line_info[2])])])
            snarl_list.append(snarl_id)
            snarl_son_dict[snarl_id] = set()

for line in openfile(snarl_json_file):
    snarl_json = json.loads(line)
    if "parent" in snarl_json.keys():
        snarl_id = '_'.join([str(x) for x in sorted([int(snarl_json['parent']['start']['node_id']), int(snarl_json['parent']['end']['node_id'])])])
        if snarl_id in snarl_son_dict:
            snarl_son_dict[snarl_id].add(snarl_json['start']['node_id'])
            snarl_son_dict[snarl_id].add(snarl_json['end']['node_id'])


node_len_dict = {}
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        node_len_dict[node] = len(line.strip().split()[2])

fo1 = open(out_snarl_bed_file, "w")
for line in openfile(gfa_file):
    if line[0] == "W":
        contig = line.strip().split()[3]
        contig_start = int(line.strip().split()[4])
        contig_end = int(line.strip().split()[5])
        path_info = line.strip().split()[6]
        path_node = re.split('<|>', path_info)[1:]
        clip_path = False
        snarl_select_list = []
        for snarl in snarl_list:
            if snarl.split("_")[0] in path_node or snarl.split("_")[1] in path_node:
                snarl_select_list.append([snarl.split("_")[0], snarl.split("_")[1]])
                clip_path = True
        if clip_path:
            path_pos = []
            pos = contig_start
            for node in path_node:
                path_pos.append(pos)
                pos += node_len_dict[node]
            for snarl in snarl_select_list:
                snarl_id = '_'.join([str(x) for x in sorted([int(snarl[0]), int(snarl[1])])])
                son_nodes = snarl_son_dict[snarl_id]
                node_index = [index for index, node in enumerate(path_node) if node == snarl[0] or node == snarl[1]]
                for i, index in enumerate(node_index):
                    if i == 0:
                        previous_index = 0
                    else:
                        previous_index = node_index[i - 1]
                    if i == len(node_index) - 1:
                        next_index = len(path_node) - 1
                    else:
                        next_index = node_index[i + 1]
                    left_son_nodes = set(path_node[previous_index + 1 : index]) & son_nodes
                    right_son_nodes = set(path_node[index + 1 : next_index]) & son_nodes
                    if len(left_son_nodes) > len(right_son_nodes):
                        candidate_nodes = left_son_nodes | {snarl[0], snarl[1]}
                        for j in list(range(previous_index, index)):
                            if path_node[j] in candidate_nodes:
                                fo1.write(contig + '\t' + str(path_pos[j]) + '\t' + str(path_pos[index] + node_len_dict[path_node[index]]) + '\n')
                                break
                    elif len(left_son_nodes) < len(right_son_nodes):
                        candidate_nodes = right_son_nodes | {snarl[0], snarl[1]}
                        for j in list(range(index + 1, next_index + 1))[::-1]:
                            if path_node[j] in candidate_nodes:
                                fo1.write(contig + '\t' + str(path_pos[index]) + '\t' + str(path_pos[j] + node_len_dict[path_node[j]]) + '\n')
                                break
                    fo1.write(contig + '\t' + str(max(path_pos[index] - boundary_size, contig_start))  + '\t' + str(min(path_pos[index] + node_len_dict[path_node[index]] + boundary_size, contig_end)) + '\n')
fo1.close()

