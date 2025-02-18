#!/usr/bin/env python3
import sys
import gzip
import re
from collections import Counter
import networkx as nx


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

def main():
    gfa_file = sys.argv[1]
    context_len = int(sys.argv[2])
    ref_num = int(sys.argv[3])
    exten_size = 10
    max_region_size = 10000
    jaccard_similarity = 0
    ref_include_sample_file = sys.argv[4]
    ref_candidate_sample_file = sys.argv[5]
    out_gfa_file = sys.argv[6]
    out_region_file = sys.argv[7]

    ref_include_sample_list = [line.strip() for line in openfile(ref_include_sample_file)]
    ref_candidate_sample_list = [line.strip() for line in openfile(ref_candidate_sample_file)]

    print("0. Load the node info.")
    node_index = 0
    node_len_dict = {}
    for line in openfile(gfa_file):
        if line[0] == "S":
            node = line.strip().split()[1]
            node_len_dict[node] = len(line.strip().split()[2])
            node_index = int(node) if int(node) > node_index else node_index

    print("1. Select reference haplotypes.")
    ######Detect the cycle nodes in the reference haplotypes
    hap_cycle_len_dict = {}
    ref_hap_list = []
    processing = 0
    for line in openfile(gfa_file):
        if line[0] == "W":
            sample = line.strip().split()[1]
            hap = line.strip().split()[2]
            hap_name = sample + "#" + hap
            if sample in ref_include_sample_list:
                if hap_name not in ref_hap_list:
                    ref_hap_list.append(hap_name)
            if sample in ref_candidate_sample_list:
                path_info = line.strip().split()[6]
                path_node = re.split('<|>', path_info)[1:]
                path_node_count = Counter(path_node)
                cycle_len = 0
                for node, count in path_node_count.items():
                    if count >= 2:
                        cycle_len += count * node_len_dict[node]
                hap_cycle_len_dict[hap_name] = hap_cycle_len_dict.get(hap_name, 0) + cycle_len
            processing += 1
            print("W line processing {}".format(processing), end="\r")
            if processing % 100 == 0:
                print("W line processing {}".format(processing), end="\r")
    sorted_hap_cycle_len_list = [[info[0], info[1]] for info in sorted(hap_cycle_len_dict.items(), key=lambda x: x[1], reverse=True)]

    ref_current_num = 0
    for info in sorted_hap_cycle_len_list:
        sample = info[0].split("#")[0]
        ref_hap_list.append(info[0])
        ref_current_num += 1
        if ref_current_num >= ref_num:
            break

    print("")
    print("Select {} reference haplotypes with the longest cycle region are selected: {}".format(len(ref_hap_list), "\t".join(ref_hap_list)))

    print("2. Detect the cycle nodes.")
    cycle_node_dict = {}
    for line in openfile(gfa_file):
        if line[0] == "W":
            sample = line.strip().split()[1]
            hap = line.strip().split()[2]
            hap_name = sample + "#" + hap
            if hap_name in ref_hap_list:
                path_info = line.strip().split()[6]
                path_node = re.split('<|>', path_info)[1:]
                path_node_count = Counter(path_node)
                for node, count in path_node_count.items():
                    if count >= 2:
                        cycle_node_dict[node] = True
            processing += 1
            if processing % 100 == 0:
                print("W line processing {}".format(processing), end="\r")
    print("{} cycle nodes are found.".format(len(cycle_node_dict)))

    print("3. Detect the cycle regions in the reference haplotypes.")
    ######Detect the cycle regions in the each reference haplotypes
    node_region_dict = {}
    for node in cycle_node_dict:
        node_region_dict[node] = []
    region_info_dict = {}
    region_index = 0
    processing = 0
    for line in openfile(gfa_file):
        if line[0] == "W":
            sample = line.strip().split()[1]
            hap = line.strip().split()[2]
            contig = line.strip().split()[3]
            hap_name = sample + "#" + hap
            if hap_name in ref_hap_list:
                path_info = line.strip().split()[6]
                path_node = re.split('<|>', path_info)[1:]
                ####detect the cycle region
                cycle_region_list = []
                pos = int(line.strip().split()[4])
                start_pos = int(line.strip().split()[4])
                start_i = 0
                end_i = 0
                non_cycle_len = 0
                is_cycle = False
                for i in range(len(path_node)):
                    node = path_node[i]
                    if not is_cycle:
                        if node in cycle_node_dict:
                            is_cycle = True
                            start_i = i
                            start_pos = pos
                            non_cycle_len = 0
                    elif is_cycle:
                        if node in cycle_node_dict:
                            non_cycle_len = 0
                        else:
                            if non_cycle_len == 0:
                                end_i = i
                            non_cycle_len += node_len_dict[node]
                            if non_cycle_len > exten_size:
                                is_cycle = False
                                node_list = path_node[start_i: end_i]
                                first_pos = start_pos
                                second_pos = start_pos
                                first_i = start_i
                                second_i = start_i
                                for j in range(len(node_list)):
                                    second_pos += node_len_dict[node_list[j]]
                                    second_i += 1
                                    if second_pos - first_pos >= max_region_size:
                                        cycle_region_list.append([first_i, second_i, first_pos, second_pos])
                                        first_i = second_i
                                        first_pos = second_pos
                                if second_pos != first_pos:
                                    cycle_region_list.append([first_i, second_i, first_pos, second_pos])
                    pos += node_len_dict[node]
                if is_cycle:
                    if non_cycle_len == 0:
                        end_i = i + 1
                    node_list = path_node[start_i: end_i]
                    first_pos = start_pos
                    second_pos = start_pos
                    first_i = start_i
                    second_i = start_i
                    for j in range(len(node_list)):
                        second_pos += node_len_dict[node_list[j]]
                        second_i += 1
                        if second_pos - first_pos >= max_region_size:
                            cycle_region_list.append([first_i, second_i, first_pos, second_pos])
                            first_i = second_i
                            first_pos = second_pos
                    if second_pos != first_pos:
                        cycle_region_list.append([first_i, second_i, first_pos, second_pos])
                
                for cycle_region in cycle_region_list:
                    first_i = cycle_region[0]
                    second_i = cycle_region[1]
                    first_pos = cycle_region[2]
                    second_pos = cycle_region[3]
                    ####Node region mapping
                    region_node = set(path_node[first_i: second_i])
                    for node in region_node:
                        if node in cycle_node_dict:
                            node_region_dict[node].append(region_index)
                    ####Region context mapping
                    context_start = max(0, first_i - context_len)
                    context_end = min(len(path_node), second_i + context_len)
                    context_node = set(path_node[context_start: context_end])
                    region_info_dict[region_index] = [hap_name, region_node, context_node, contig, [first_pos, second_pos]]
                    region_index += 1

            processing += 1
            if processing % 100 == 0:
                print("W line processing {}".format(processing), end="\r")
    print("")
    print("{} cycle regions are found.".format(region_index))


    print("4. Construct the orthology graph.")
    region_match_dict = {}
    processing = 0
    ######Construct the orthology graph among different reference haplotypes
    for query_region, query_info in region_info_dict.items():
        query_hap_name = query_info[0]
        query_region_node = query_info[1]
        query_context_node = query_info[2]
        query_pos = query_info[4]
        query_len = query_pos[1] - query_pos[0]
        target_region_candidate_list = []
        ####Select the target regions based on overlapped node
        for node in query_region_node:
            if node in cycle_node_dict:
                for target_region in node_region_dict[node]:
                    target_region_candidate_list.append(target_region)
        target_region_candidate_set = set(target_region_candidate_list)
        ####Calculate the Jaccord similarity between different regions
        hap_jaccard_sim_dict = {}
        query_context_node_len = len(query_context_node)
        for target_region in target_region_candidate_set:
            target_info = region_info_dict[target_region]
            target_hap_name = target_info[0]
            if target_hap_name == query_hap_name:
                continue
            target_context_node = target_info[2]
            target_contig = target_info[3]
            target_pos = target_info[4]
            target_context_node_len = len(target_context_node)
            jaccard_context_node_isec_len = len(query_context_node & target_context_node)
            jaccard_sim = jaccard_context_node_isec_len/(query_context_node_len + target_context_node_len - jaccard_context_node_isec_len)
            if jaccard_sim > jaccard_similarity:
                if target_hap_name not in hap_jaccard_sim_dict:
                    hap_jaccard_sim_dict[target_hap_name] = [[target_region, jaccard_sim, target_contig, target_pos]]
                else:
                    hap_jaccard_sim_dict[target_hap_name].append([target_region, jaccard_sim, target_contig, target_pos])
        
        region_match_dict[query_region] = []
        for target_hap_name, hap_jaccard_sim_list in hap_jaccard_sim_dict.items():
            if len(hap_jaccard_sim_list) != 0:
                sorted_hap_jaccard_sim_list = sorted(hap_jaccard_sim_list, key = lambda x: x[1], reverse = True)
                ####find the best pos
                best_region = sorted_hap_jaccard_sim_list[0][0]
                best_contig = sorted_hap_jaccard_sim_list[0][2]
                best_pos = sorted_hap_jaccard_sim_list[0][3]
                isec_len = len(region_info_dict[best_region][1] & query_region_node)
                region_match_dict[query_region].append(best_region)
                ####find the next best pos with in the best target pos plus query region length
                for i in range(1, len(sorted_hap_jaccard_sim_list)):
                    pos_dist = max(sorted_hap_jaccard_sim_list[i][3][0] - best_pos[1], best_pos[0] - sorted_hap_jaccard_sim_list[i][3][1])
                    if pos_dist <= query_len - isec_len and best_contig == sorted_hap_jaccard_sim_list[i][2]:
                        region_match_dict[query_region].append(sorted_hap_jaccard_sim_list[i][0])
                    else:
                        break

        processing += 1
        if processing % 100 == 0:
            print("Region processing {}".format(processing), end="\r")
    print("")
    ####Create the reciprocally matching region graph
    region_G = nx.Graph()
    for region in region_info_dict:
        region_G.add_node(region)
    for query_region, target_region_list in region_match_dict.items():
        for target_region in target_region_list:
            if query_region in region_match_dict[target_region]:
                region_G.add_edge(query_region, target_region)
    ####Detect the connect components
    subgraph_index = 0
    region_subgraph_dict = {}
    subgraph_node_sub_dict = {}
    for region_sub_G in nx.connected_components(region_G):
        node_sub_dict = {}
        node_sub_set = set()
        for region in region_sub_G:
            region_subgraph_dict[region] = subgraph_index
            node_sub_set = node_sub_set | region_info_dict[region][1]
        for node in node_sub_set:
            if node in cycle_node_dict:
                node_index += 1
                node_sub_dict[node] = str(node_index)
        subgraph_node_sub_dict[subgraph_index] = node_sub_dict
        subgraph_index += 1
    print("{} connect connect components are found.".format(subgraph_index))
    fo = open(out_region_file, "w")
    for region, info in region_info_dict.items():
        fo.write('\t'.join([info[0], info[3], str(info[4][0]), str(info[4][1]), str(region), str(region_subgraph_dict[region])] + [','.join(list(info[1]))] ) + "\n")
    fo.close()
    print("{} connect connect components are found.".format(subgraph_index))

    print("5. Load the cycle node and edge info")
    ######Load the cycle node and edge info
    fo = open(out_gfa_file, "w")
    cycle_node_info_dict = {}
    cycle_edge_info_dict = {}
    for line in openfile(gfa_file):
        if line[0] == "S":
            node_info = line.strip().split()[1:]
            node = node_info[0]
            if node in cycle_node_dict:
                cycle_node_info_dict[node] = node_info
            fo.write(line)
        elif line[0] == "L":
            edge_info = line.strip().split()[1:]
            edge = edge_info[0] + edge_info[1] + edge_info[2] + edge_info[3]
            node1 = edge_info[0]
            node2 = edge_info[2]
            if node1 in cycle_node_dict or node2 in cycle_node_dict:
                cycle_edge_info_dict[edge] = edge_info
            fo.write(line)


    print("6. Cycle node, edge and path substitution")
    ######Cycle node, edge and path substitution
    processing = 0
    sub_cycle_node_info_dict = {}
    sub_cycle_edge_info_dict = {}
    for line in openfile(gfa_file):
        if line[0] == "W":
            sample = line.strip().split()[1]
            hap = line.strip().split()[2]
            hap_name = sample + "#" + hap
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            path_size = len(path_node)
            path_cycle_node_dict = {}
            sub_path_node = path_node.copy()
            cycle_region_list = []
            start_i = 0
            end_i = 0
            non_cycle_len = 0
            is_cycle = False
            for i in range(path_size):
                node = path_node[i]
                if not is_cycle:
                    if node in cycle_node_dict:
                        path_cycle_node_dict[node] = True
                        is_cycle = True
                        start_i = i
                        non_cycle_len = 0
                elif is_cycle:
                    if node in cycle_node_dict:
                        path_cycle_node_dict[node] = True
                        non_cycle_len = 0
                    else:
                        if non_cycle_len == 0:
                            end_i = i
                        non_cycle_len += node_len_dict[node]
                        if non_cycle_len > exten_size:
                            is_cycle = False
                            node_list = path_node[start_i: end_i]
                            first_pos = 0
                            second_pos = 0
                            first_i = start_i
                            second_i = start_i
                            for j in range(len(node_list)):
                                second_pos += node_len_dict[node_list[j]]
                                second_i += 1
                                if second_pos - first_pos >= max_region_size:
                                    cycle_region_list.append([first_i, second_i])
                                    first_i = second_i
                                    first_pos = second_pos
                            if second_pos != first_pos:
                                cycle_region_list.append([first_i, second_i])
            if is_cycle:
                if non_cycle_len == 0:
                    end_i = i + 1
                node_list = path_node[start_i: end_i]
                first_pos = 0
                second_pos = 0
                first_i = start_i
                second_i = start_i
                for j in range(len(node_list)):
                    second_pos += node_len_dict[node_list[j]]
                    second_i += 1
                    if second_pos - first_pos >= max_region_size:
                        cycle_region_list.append([first_i, second_i])
                        first_i = second_i
                        first_pos = second_pos
                if second_pos != first_pos:
                    cycle_region_list.append([first_i, second_i])

            for cycle_region in cycle_region_list:
                first_i = cycle_region[0]
                second_i = cycle_region[1]
                ####Node region mapping
                region_node = set(path_node[first_i: second_i])
                context_start = max(0, first_i - context_len)
                context_end = min(path_size, second_i + context_len)
                context_node = set(path_node[context_start: context_end])
                match_region = ""
                match_jaccard_sim = 0
                target_region_candidate_list = []
                ####Search the target region
                for node in region_node:
                    if node in path_cycle_node_dict:
                        for target_region in node_region_dict[node]:
                            target_region_candidate_list.append(target_region)
                target_region_candidate_set = set(target_region_candidate_list)
                ####Calculate the best matching region
                context_node_len = len(context_node)
                for target_region in target_region_candidate_set:
                    target_info = region_info_dict[target_region]
                    target_context_node = target_info[2]
                    target_context_node_len = len(target_context_node)
                    jaccard_context_node_isec_len = len(context_node & target_context_node)
                    jaccard_sim = jaccard_context_node_isec_len/(context_node_len + target_context_node_len - jaccard_context_node_isec_len)
                    if jaccard_sim > match_jaccard_sim:
                        match_region = target_region
                        match_jaccard_sim = jaccard_sim
                if match_region != "":
                    subgraph = region_subgraph_dict[match_region]
                    node_sub_dict = subgraph_node_sub_dict[subgraph]
                    ####Subgraph substitution
                    for node in region_node:
                        if node not in node_sub_dict:
                            if node in path_cycle_node_dict:
                                node_index += 1
                                node_sub_dict[node] = str(node_index)
                    
                    ##path substitution
                    for i in range(first_i, second_i):
                        node = path_node[i]
                        if node in path_cycle_node_dict:
                            sub_path_node[i] = node_sub_dict[node]
                    
            for i in range(path_size):
                ##node info
                if path_node[i] != sub_path_node[i]:
                    node = path_node[i]
                    node_info = cycle_node_info_dict[node].copy()
                    node_info[0] = sub_path_node[i]
                    sub_cycle_node_info_dict[node_info[0]] = node_info
                ##edge info
                if i != path_size - 1:
                    if path_node[i] != sub_path_node[i] or path_node[i + 1] != sub_path_node[i + 1]:
                        edge = path_node[i] + path_direct[i].replace(">","+").replace("<","-") + path_node[i + 1] + path_direct[i + 1].replace(">","+").replace("<","-")
                        rev_edge = rev_edge_info(edge)
                        if edge in cycle_edge_info_dict:
                            edge_info = cycle_edge_info_dict[edge].copy()
                            edge_info[0] = sub_path_node[i]
                            edge_info[2] = sub_path_node[i + 1]
                        elif rev_edge in cycle_edge_info_dict:
                            edge_info = cycle_edge_info_dict[rev_edge].copy()
                            edge_info[0] = sub_path_node[i + 1]
                            edge_info[2] = sub_path_node[i]
                        else:
                            print("{} is not found!".format(''.join(edge_info[0:4])))
                        sub_cycle_edge_info_dict[''.join(edge_info[0:4])] = edge_info
          
            ##write path
            write_path = ''.join([path_direct[j] + sub_path_node[j] for j in range(path_size)])
            fo.write('\t'.join(["W"] + line.strip().split()[1:6]) + "\t" + write_path + '\n')
            
            processing += 1
            if processing % 100 == 0:
                print("W line processing {}".format(processing), end="\r")
    print("")
    ##write node
    for node_info in sub_cycle_node_info_dict.values():
        fo.write('\t'.join(["S"] + node_info) + '\n')
    ##write edge
    for edge_info in sub_cycle_edge_info_dict.values():
        fo.write('\t'.join(["L"] + edge_info) + '\n')
    fo.close()

if __name__ == "__main__":
    main()
