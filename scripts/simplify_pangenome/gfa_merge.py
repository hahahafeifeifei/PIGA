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


gfa_list = sys.argv[1]
out_gfa_file = sys.argv[2]
prefix = sys.argv[3]
ref_sample = sys.argv[4]
threads = int(sys.argv[5])

ref_pos_list = []
for gfa_file in openfile(gfa_list):
    gfa_file = gfa_file.strip()
    for line in openfile(gfa_file):
        line_info = line.strip().split()
        if line[0] == "W" and line.split()[1] == ref_sample:
            ref_pos_list.append([gfa_file, line_info[1], int(line_info[4])])

ref_sorted_pos_list = sorted(ref_pos_list, key=lambda x: x[2])
gfa_index_dict = {}
for i, item in enumerate(ref_sorted_pos_list):
    key = item[0]
    if key not in gfa_index_dict:
        gfa_index_dict[key] = i


fo1 = open(prefix + '.S_L', "w")
fo2 = open(prefix + '.W_subgraph', "w")
for gfa_file in openfile(gfa_list):
    gfa_file = gfa_file.strip()
    for line in openfile(gfa_file):
        if line[0] == "S" or line[0] == "L":
            fo1.write(line)
        elif line[0] == "W":
            fo2.write("\t".join(line.strip().split() + [str(gfa_index_dict[gfa_file])]) + "\n")
fo1.close()
fo2.close()

subprocess.run("sort -k1,1 -k2,2 -k3,3n -k4,4 -k5,5n -T . --buffer-size=60G --parallel={} {} > {}".format(threads, prefix + ".W_subgraph", prefix + ".sort.W_subgraph"), shell = True, stdout = None, stderr = None)
subprocess.run("rm {}".format(prefix + ".W_subgraph"), shell = True, stdout = None, stderr = None)
pre_contig = "null"
add_edge_dict = {}
fo = open(prefix + ".connect.W", "w")
for line in openfile(prefix + ".sort.W_subgraph"):
    line_info = line.strip().split()
    sample = line_info[1]
    hap = line_info[2]
    contig = line_info[3]
    start = line_info[4]
    end = line_info[5]
    path = line_info[6]
    subgraph = int(line_info[7])
    if pre_contig == "null":
        pre_sample = sample
        pre_hap = hap
        pre_contig = contig
        pre_path = path
        pre_start = start
        pre_end = end
        pre_subgraph = subgraph
    elif pre_sample == sample and contig == pre_contig and start == pre_end and abs(subgraph-pre_subgraph)<=1:
        pre_path_end_node = re.split('<|>', pre_path)[1:][-1]
        pre_path_end_direct = re.split('\d+', pre_path)[:-1][-1].replace(">","+").replace("<","-")
        path_start_node = re.split('<|>', path)[1:][0]
        path_start_direct = re.split('\d+', path)[:-1][0].replace(">","+").replace("<","-")
        add_edge_dict['\t'.join([pre_path_end_node, pre_path_end_direct, path_start_node, path_start_direct])] = True
        pre_end = end
        pre_subgraph = subgraph
        pre_path += path
    else:
        fo.write('\t'.join(["W", pre_sample, pre_hap, pre_contig, pre_start, pre_end, pre_path]) + '\n')
        pre_sample = sample
        pre_hap = hap
        pre_contig = contig
        pre_path = path
        pre_start = start
        pre_end = end
        pre_subgraph = subgraph
fo.write('\t'.join(["W", pre_sample, pre_hap, pre_contig, pre_start, pre_end, pre_path]) + '\n')
fo.close()
subprocess.run("rm {}".format(prefix + ".sort.W_subgraph"), shell = True, stdout = None, stderr = None)

fo = open(out_gfa_file, "w")
for line in openfile(prefix + ".S_L"):
    if line[0] == "S":
        fo.write(line)

for line in openfile(prefix + ".S_L"):
    if line[0] == "L":
        fo.write(line)
for edge in add_edge_dict:
    fo.write("\t".join(["L", edge, "0M"]) + '\n')

for line in openfile(prefix + ".connect.W"):
    fo.write(line)
fo.close()
subprocess.run("rm {} {}".format(prefix + ".S_L", prefix + ".connect.W"), shell = True, stdout = None, stderr = None)
