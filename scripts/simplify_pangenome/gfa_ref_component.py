import sys
import gzip
import re
import networkx as nx

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

gfa_file = sys.argv[1]
out_gfa_file = sys.argv[2]
ref_sample_list = sys.argv[3].split(",")

G = nx.Graph()
for line in openfile(gfa_file):
    if line[0] == "S":
        G.add_node(int(line.strip().split()[1]))
    elif line[0] == "L":
        G.add_edge(int(line.strip().split()[1]), int(line.strip().split()[3]))

node_component_dict = {}
for component_index, component_G in enumerate(nx.connected_components(G)):
    for node in component_G:
        node_component_dict[node] = component_index

ref_component_list = []
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.strip().split()[1]
        if sample in ref_sample_list:
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            ref_component_list.append(node_component_dict[int(path_node[0])])

fo = open(out_gfa_file, "w")
for line in openfile(gfa_file):
    if line[0] == "S":
        component = node_component_dict[int(line.strip().split()[1])]
        if component in ref_component_list:
            fo.write(line)
    elif line[0] == "L":
        component1 = node_component_dict[int(line.strip().split()[1])]
        component2 = node_component_dict[int(line.strip().split()[3])]
        if component1 in ref_component_list and component2 in ref_component_list:
            fo.write(line)
    elif line[0] == "W":
        path_info = line.strip().split()[6]
        path_node = re.split('<|>', path_info)[1:]
        component = node_component_dict[int(path_node[0])]
        if component in ref_component_list:
            fo.write(line)
    else:
        fo.write(line)
fo.close()