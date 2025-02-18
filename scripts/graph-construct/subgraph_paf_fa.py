import sys
import os
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def takeFirst(elem):
    return elem[0]

def extract_subgraph_info(aln_info_list, node_subgraph_dict):
    aln_info_list.sort(key = takeFirst)
    aln_path_list = []
    aln_gl_list = []
    for aln_info in aln_info_list:
        aln_path_list += aln_info[1].replace("<",">").split(">")[1:]
        aln_gl_list += [aln_info[2]] * len(aln_info[1].replace("<",">").split(">")[1:])
    subgraph_info_list = [[], [], [], [], []]
    for i in range(len(aln_path_list)):
        node = aln_path_list[i]
        gl = aln_gl_list[i]
        subgraph = node_subgraph_dict[node]
        if subgraph_info_list == [[]] * 5:
            subgraph_info_list[0].append(subgraph)
            subgraph_info_list[1].append(node)
        elif subgraph != subgraph_info_list[0][-1]:
            subgraph_info_list[2].append(last_node)
            subgraph_info_list[3].append(last_gl)
            subgraph_info_list[0].append(subgraph)
            subgraph_info_list[1].append(node)
        last_node = node
        last_gl = gl
    subgraph_info_list[2].append(last_node)
    subgraph_info_list[3].append(last_gl)
    subgraph_info_list[4] = [""] * len(subgraph_info_list[0])
    return subgraph_info_list

def filter_subgraph_info(contig_subgraph_info):
    filtered_contig_subgraph_info = [[], [], [], [], []]
    for i in range(len(contig_subgraph_info[0])):
        if contig_subgraph_info[4][i] != "":
            subgraph = contig_subgraph_info[0][i]
            source_node = contig_subgraph_info[1][i]
            if filtered_contig_subgraph_info == [[]] * 5:
                filtered_contig_subgraph_info[0].append(subgraph)
                filtered_contig_subgraph_info[1].append(source_node)
            elif subgraph != filtered_contig_subgraph_info[0][-1]:
                filtered_contig_subgraph_info[2].append(last_sink_node)
                filtered_contig_subgraph_info[3].append(last_gl)
                filtered_contig_subgraph_info[4].append(last_end)
                filtered_contig_subgraph_info[0].append(subgraph)
                filtered_contig_subgraph_info[1].append(source_node)
            last_sink_node = contig_subgraph_info[2][i]
            last_gl = contig_subgraph_info[3][i]
            last_end = contig_subgraph_info[4][i]
    if len(filtered_contig_subgraph_info[0]) != 0:
        filtered_contig_subgraph_info[2].append(last_sink_node)
        filtered_contig_subgraph_info[3].append(last_gl)
        filtered_contig_subgraph_info[4].append(last_end)
    return filtered_contig_subgraph_info


prefix = sys.argv[1]
paf_file = sys.argv[2]
gaf_file = sys.argv[3]
outdir = sys.argv[4]

####divide node into subgraph
if "/" not in prefix:
    prefix_str = prefix
    prefix_dir = "./"
else:
    prefix_dir = "/".join(prefix.split("/")[:-1])
    prefix_str = prefix.split("/")[-1]
subgraph_file = []
for file in os.listdir(prefix_dir):
    if prefix_str in file and "gfa" in file:
        subgraph_file.append(file)

node_subgraph_dict = {}
for i in range(len(subgraph_file)):
    for line in open(prefix_dir + "/" + prefix_str + "_" + str(i) + '.gfa'):
        if line[0] == "S":
            node_id = "s" + line.strip().split()[1]
            node_subgraph_dict[node_id] = i

####search the subgraph path of each contig
contig_subgraph_info_dict = {}
last_contig = ""
aln_info_list = []
for line in openfile(gaf_file):
    contig = line.split()[0]
    if contig != last_contig and last_contig != "":
        contig_subgraph_info_dict[last_contig] = extract_subgraph_info(aln_info_list, node_subgraph_dict)
        aln_info_list = []
    aln_start = int(line.split()[2])
    aln_path = line.split()[5]
    aln_gl = line.split()[10]
    aln_info_list.append([aln_start, aln_path, aln_gl])
    last_contig = contig
contig_subgraph_info_dict[last_contig] = extract_subgraph_info(aln_info_list, node_subgraph_dict)

####find the breakpoint position in paf file
for line in openfile(paf_file):
    contig = line.split()[0]
    node = line.split()[5].split("id=_MINIGRAPH_|")[1]
    contig_subgraph_info = contig_subgraph_info_dict[contig]
    if len(contig_subgraph_info[0]) == 1:
        contig_end = int(line.split()[1])
        contig_subgraph_info[4][0] = contig_end
    elif node in contig_subgraph_info[2]:
        aln_lg = line.split()[14].split("gl:i:")[1]
        for i in range(len(contig_subgraph_info[0])):
            if contig_subgraph_info[2][i] == node and contig_subgraph_info[3][i] == aln_lg:
                aln_end = int(line.split()[3])
                contig_subgraph_info[4][i] = aln_end
        contig_end = int(line.split()[1])
        contig_subgraph_info[4][i] = contig_end
 
####delete the subgraph with no paf record (bed alignment)
delete_contig_list = []
for contig, contig_subgraph_info in contig_subgraph_info_dict.items():
    ##delete the subgraph with end pos <= previous end pos
    for i in range(len(contig_subgraph_info[0])):
        last_aln_end = 0
        for j in range(0, i):
            if contig_subgraph_info[4][j] != "":
                last_aln_end = contig_subgraph_info[4][j]
        if contig_subgraph_info[4][i] != "" and contig_subgraph_info[4][i] <= last_aln_end:
            contig_subgraph_info[4][i] = ""

    if "" in contig_subgraph_info[4]:
        filtered_contig_subgraph_info = filter_subgraph_info(contig_subgraph_info)
        if len(filtered_contig_subgraph_info[0]) == 0:
            delete_contig_list.append(contig)
        else:
            contig_subgraph_info_dict[contig] = filtered_contig_subgraph_info
for contig in delete_contig_list:
    del(contig_subgraph_info_dict[contig])

####generate the contig coordinate bed for each subgraph
fo_list = []
for i in range(len(subgraph_file)):
    fo_list.append(open(os.path.join(outdir, prefix_str + '_' + str(i) + '.bed'), "w"))
for contig, contig_subgraph_info in contig_subgraph_info_dict.items():
    sample = contig.split("id=")[1].split("|")[0]
    contig_name = contig.split("id=")[1].split("|")[1]
    for i in range(len(contig_subgraph_info[0])):
        subgraph = contig_subgraph_info[0][i]
        if i == 0:
            fo_list[subgraph].write(sample + '\t' + contig_name + ":" + str(1) + "-" + str(contig_subgraph_info[4][i]) + '\n')
        else:
            fo_list[subgraph].write(sample + '\t' + contig_name + ":" + str(contig_subgraph_info[4][i - 1] + 1) + "-" + str(contig_subgraph_info[4][i]) + '\n')
for i in range(len(subgraph_file)):
    fo_list[i].close()

####generate the paf file for each subgraph
fo_list = []
for i in range(len(subgraph_file)):
    fo_list.append(open(os.path.join(outdir, prefix_str + '_' + str(i) + '.paf'), "w"))
for line in openfile(paf_file):
    line_info = line.strip().split()
    contig = line_info[0]
    aln_subgraph = node_subgraph_dict[line_info[5].split("|")[1]]
    aln_start = int(line_info[2]) + 1
    aln_end = int(line_info[3])
    if contig not in contig_subgraph_info_dict.keys():
        continue
    contig_subgraph_info = contig_subgraph_info_dict[contig]
    for i in range(len(contig_subgraph_info[0])):
        subgraph = contig_subgraph_info[0][i]
        if i == 0:
            contig_start = 1
            contig_end = contig_subgraph_info[4][i]
        else:
            contig_start = contig_subgraph_info[4][i - 1] + 1
            contig_end = contig_subgraph_info[4][i]
        if contig_start <= aln_start and contig_end >= aln_end and aln_subgraph == subgraph:
            contig_split = contig + ":" + str(contig_start) + "-" + str(contig_end)
            aln_split_start = aln_start - contig_start + 1
            aln_split_end = aln_end - contig_start + 1
            line_info[0] = contig_split
            line_info[2] = str(aln_split_start - 1)
            line_info[3] = str(aln_split_end)
            fo_list[subgraph].write("\t".join(line_info) + '\n')
            break

for i in range(len(subgraph_file)):
    fo_list[i].close()

