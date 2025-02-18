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

def interval_overlap(start1, end1, start2, end2):
    max_start = start1 if start1 >= start2 else start2
    min_end = end1 if end1 <= end2 else end2
    overlap_len = min_end - max_start + 1
    return True if overlap_len > 0 else False

def get_edge_seq(seq, direct):
    if direct == ">":
        return seq
    elif direct == "<":
        rev_dict = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
        rev_seq = ''.join([rev_dict[base] for base in seq[::-1]])
        return rev_seq


gfa_file = sys.argv[1]
bed_file = sys.argv[2]
ref_sample_list = sys.argv[3].split(",")
output_file = sys.argv[4]

clip_region_dict = {}
for line in openfile(bed_file):
    contig = line.strip().split()[0]
    start = int(line.strip().split()[1])
    end = int(line.strip().split()[2])
    if contig not in clip_region_dict:
        clip_region_dict[contig] = [[start, end]]
    else:
        clip_region_dict[contig].append([start, end])

node_seq_dict = {}
node_index = 0
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        node_index = node_index if node_index >= int(node) else int(node)
        node_seq_dict[node] = line.strip().split()[2]

fo = open(output_file, "w")
for line in openfile(gfa_file):
    if line[0] == "S":
        fo.write(line)
    elif line[0] == "L":
        fo.write(line)
    elif line[0] == "W":
        contig = line.strip().split()[3]
        if contig in clip_region_dict:
            sample = line.strip().split()[1]
            hap = line.strip().split()[2]
            contig_start = int(line.strip().split()[4])
            contig_end = int(line.strip().split()[5])
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_i = 0

            if sample in ref_sample_list:
                write_path = ""
                pos = contig_start
                for i in range(len(path_node)):
                    node = path_node[i]
                    direct = path_direct[i]
                    node_end = pos + len(node_seq_dict[node])
                    is_clip = False
                    for mask in clip_region_dict[contig]:
                        mask_start = mask[0]
                        mask_end = mask[1]
                        if interval_overlap(pos + 1, node_end, mask_start + 1, mask_end):
                            is_clip = True
                            break
                    
                    if not is_clip:
                        if pre_i != i:
                            seq = ''.join([get_edge_seq(node_seq_dict[path_node[j]], path_direct[j]) for j in range(pre_i, i)])
                            node_index += 1
                            write_path += ">" + str(node_index)
                            fo.write('\t'.join(["S", str(node_index), seq]) + '\n')
                            if pre_i > 0:
                                fo.write('\t'.join(["L", path_node[pre_i - 1], path_direct[pre_i - 1].replace(">","+").replace("<","-"), str(node_index), "+", "0M"]) + '\n')
                            fo.write('\t'.join(["L", str(node_index), "+", path_node[i], path_direct[i].replace(">","+").replace("<","-"), "0M"]) + '\n')
                        write_path += path_direct[i] + path_node[i]
                        pre_i = i + 1
                    pos = node_end

                if pre_i != i + 1:
                    seq = ''.join([get_edge_seq(node_seq_dict[path_node[j]], path_direct[j]) for j in range(pre_i, i + 1)])
                    node_index += 1
                    write_path += ">" + str(node_index)
                    fo.write('\t'.join(["S", str(node_index), seq]) + '\n')
                    if pre_i > 0:
                        fo.write('\t'.join(["L", path_node[pre_i - 1], path_direct[pre_i - 1].replace(">","+").replace("<","-"), str(node_index), "+", "0M"]) + '\n')
                fo.write('\t'.join(["W", sample, hap, contig, str(contig_start), str(contig_end), write_path]) + '\n')

            else:
                pre_pos = contig_start
                pos = contig_start
                for i in range(len(path_node)):
                    node = path_node[i]
                    direct = path_direct[i]
                    node_info = node + direct.replace(">","+").replace("<","-")
                    node_end = pos + len(node_seq_dict[node])
                    for mask in clip_region_dict[contig]:
                        mask_start = mask[0]
                        mask_end = mask[1]
                        if interval_overlap(pos + 1, node_end, mask_start + 1, mask_end):
                            if pre_i != i:
                                write_path = ''.join([path_direct[j] + path_node[j] for j in range(pre_i, i)])
                                fo.write('\t'.join(["W", sample, hap, contig, str(pre_pos), str(pos), write_path]) + '\n')
                            pre_i = i + 1
                            pre_pos = node_end
                            break
                    pos = node_end

                if pre_i != i + 1:
                    write_path = ''.join([path_direct[j] + path_node[j]for j in range(pre_i, i + 1)])
                    fo.write('\t'.join(["W", sample, hap, contig, str(pre_pos), str(pos), write_path]) + '\n')
        else:
            fo.write(line)
fo.close()
