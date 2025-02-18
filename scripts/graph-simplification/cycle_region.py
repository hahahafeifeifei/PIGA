
import sys
import gzip
import re
from collections import Counter

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")
def main():
    gfa_file = sys.argv[1]
    out_region_file = sys.argv[2]

    print("0. Load the node info.")
    node_index = 0
    node_len_dict = {}
    for line in openfile(gfa_file):
        if line[0] == "S":
            node = line.strip().split()[1]
            node_len_dict[node] = len(line.strip().split()[2])
            node_index = int(node) if int(node) > node_index else node_index

    fo = open(out_region_file, "w")
    print("1. Detect the cycle nodes and select reference haplotypes.")
    ######Detect the cycle nodes in the reference haplotypes
    processing = 0
    for line in openfile(gfa_file):
        if line[0] == "W":
            sample = line.strip().split()[1]
            hap = line.strip().split()[2]
            contig = line.strip().split()[3]
            hap_name = sample + "#" + hap
            cycle_dict = {}
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_node_count = Counter(path_node)
            for node, count in path_node_count.items():
                if count >= 2:
                    cycle_dict[node] = True

            pos = int(line.strip().split()[4])
            start_pos = int(line.strip().split()[4])
            start_i = 0
            end_i = 0
            is_cycle = False
            for i in range(len(path_node)):
                node = path_node[i]
                if not is_cycle:
                    if node in cycle_dict:
                        is_cycle = True
                        start_i = i
                        start_pos = pos
                elif is_cycle:
                    if node not in cycle_dict:
                        is_cycle = False
                        end_pos = pos
                        end_i = i
                        fo.write(hap_name + "\t" + contig + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + ','.join(path_node[start_i: end_i]) + "\n")
                pos += node_len_dict[node]
            if is_cycle:
                end_i = i + 1
                end_pos = pos
                fo.write(hap_name + "\t" + contig + "\t" + str(start_pos) + "\t" + str(end_pos) + "\t" + ','.join(path_node[start_i: end_i]) + "\n")

            processing += 1
            if processing % 100 == 0:
                print("W line processing {}".format(processing), end="\r")
    fo.close()

if __name__ == "__main__":
    main()
