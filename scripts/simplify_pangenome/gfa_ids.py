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

def process_w_line(line_info):
    global node_dict
    path_info = line_info[6]
    path_node = re.split('<|>', path_info)[1:]
    path_direct = re.split('\d+', path_info)[:-1]
    path_node_index = [node_dict[node] for node in path_node]
    path_info_new = ''.join([path_direct[i] + path_node_index[i] for i in range(len(path_node_index))])
    line_info[6] = path_info_new
    return "\t".join(line_info) + '\n'


subgraph_index = 0
gfa_file = sys.argv[1]
out_gfa_file = sys.argv[2]
node_index = int(sys.argv[3])
fo = open(out_gfa_file, "w")
node_dict = {}

for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        node_index += 1
        node_dict[node] = str(node_index)

for line in openfile(gfa_file):
    line_info = line.strip().split()
    if line[0] == "S":
        node = line_info[1]
        line_info[1] = node_dict[node]
        fo.write("\t".join(line_info) + '\n')

    elif line[0] == "L":
        node1 = line_info[1]
        line_info[1] = node_dict[node1]
        node2 = line_info[3]
        line_info[3] = node_dict[node2]
        fo.write("\t".join(line_info) + '\n')

    elif line[0] == "W" and line_info[1]!="_MINIGRAPH_":
        result = process_w_line(line_info)
        fo.write(result)

fo.close()
