import sys
import gzip
import argparse
import bdsg
from bdsg.bdsg import PackedGraph
from pysam import VariantFile

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


def parse_args():
    parser = argparse.ArgumentParser(description='snv projection for graph', usage='snv_projection.py [option]')
    parser.add_argument('-g','--input-graph', metavar='<file>', type=str, default='None', help='Iuput graph file (VG PackedGraph)', required=True)
    parser.add_argument('-v','--input-vcf', metavar='<file>', type=str, default='None', help='Iuput vcf file', required=True)
    parser.add_argument('-r','--ref-sample', metavar='<str>', type=str, default='None', help='Input reference sample of vcf file', required=True)
    parser.add_argument('-s','--sample-suffix', metavar='<str>', type=str, default='.variant', help='The suffix of sample name in the path line')
    parser.add_argument('-p','--contig-prefix', metavar='<str>', type=str, default='variant_chain_', help='The prefix of contig name in the path line')
    parser.add_argument('-o','--output-graph', metavar='<file>', type=str, default='None', help='Output graph file (VG PackedGraph)', required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    vcf_file = VariantFile(args.input_vcf)
    sample_list = list(vcf_file.header.samples)
    g = PackedGraph()
    g.deserialize(args.input_graph)

    ref_path_list = []
    g.for_each_path_handle(lambda p: ref_path_list.append(g.get_path_name(p)) or True if g.get_path_name(p).split("#")[0] == args.ref_sample else True)

    chain_i = 0
    node_split_pos_dict = {}
    for ref_path in ref_path_list:
        path_contig = ref_path.split("#")[2]
        if path_contig not in list(vcf_file.header.contigs):
            continue
        path_range = ref_path.split("#")[3]
        path_handle = g.get_path_handle(ref_path)
        ####Select the reference path region
        if "[" in path_range and "]" in path_range and "-" in path_range:
            path_start = int(path_range[path_range.index("[") + 1 : path_range.index("]")].split("-")[0])
        else:
            path_start = 0
        path_end = path_start

        path_step = g.path_begin(path_handle)
        while(g.has_next_step(path_step)):
            path_end += g.get_length(g.get_handle_of_step(path_step))
            path_step = g.get_next_step(path_step)
        path_end += g.get_length(g.get_handle_of_step(path_step))

        
        ####Add the split pos and variant list
        split_pos_set = set()
        variant_list = []
        for variant_record in vcf_file.fetch(path_contig, path_start, path_end):
            variant_start = variant_record.pos - 1
            variant_end = variant_record.pos + len(variant_record.ref) - 1
            split_pos_set.add(variant_start)
            split_pos_set.add(variant_end)
            variant_list.append([variant_start, variant_end])
        split_pos_list = sorted(list(split_pos_set))
        print("{} variants located in {}.".format(len(variant_list), path_contig + ":" + str(path_start + 1) + "-" + str(path_end)))
        
        
        ####Detect the split pos of nodes
        path_step = g.path_begin(path_handle)
        node_start = path_start
        node_end = path_start
        split_pos = -1
        split_i = 0
        split_n = len(split_pos_list)
        while g.has_next_step(path_step):
            node_id = g.get_id(g.get_handle_of_step(path_step))
            node_start = node_end
            node_end = node_start + g.get_length(g.get_handle_of_step(path_step))
            node_split_pos = set()
            while split_pos <= node_end:
                if split_pos > node_start and split_pos < node_end:
                    if not g.get_is_reverse(g.get_handle_of_step(path_step)):
                        node_split_pos.add(split_pos - node_start)
                    else:
                        node_split_pos.add(node_end - split_pos)
                if split_i != split_n:
                    split_pos = split_pos_list[split_i]
                    split_i += 1
                else:
                    break
            if node_split_pos != set():
                if node_id not in node_split_pos_dict:
                    node_split_pos_dict[node_id] = node_split_pos
                else:
                    node_split_pos_dict[node_id] = node_split_pos_dict[node_id] | node_split_pos
            path_step = g.get_next_step(path_step)
        node_id = g.get_id(g.get_handle_of_step(path_step))
        node_start = path_start
        node_end = node_start + g.get_length(g.get_handle_of_step(path_step))
        node_split_pos = set()
        while split_pos <= node_end:
            if split_pos > node_start and split_pos < node_end:
                if not g.get_is_reverse(g.get_handle_of_step(path_step)):
                    node_split_pos.add(split_pos - node_start)
                else:
                    node_split_pos.add(node_end - split_pos)
            if split_i != split_n:
                split_pos = split_pos_list[split_i]
                split_i += 1
            else:
                break
        if node_split_pos != set():
            if node_id not in node_split_pos_dict:
                node_split_pos_dict[node_id] = node_split_pos
            else:
                node_split_pos_dict[node_id] = node_split_pos_dict[node_id] | node_split_pos
    

    ####Split the variant nodes
    split_node_map_dict = {}
    for node_id, node_split_pos in node_split_pos_dict.items():
        split_handles = g.divide_handle(g.get_handle(node_id), bdsg.std.vector_unsigned_long(sorted(list(node_split_pos))))
        split_node_map_dict[node_id] = [g.get_id(handle) for handle in split_handles]
    print("{} variant nodes are split.".format(len(node_split_pos_dict)))


    ####Patch the path with the reverse split node at path start or forward split node at path end
    path_list = []
    g.for_each_path_handle(lambda p: path_list.append(g.get_path_name(p)) or True)
    patch_i = 0
    for path in path_list:
        path_handle = g.get_path_handle(path)
        start_handle = g.path_begin(path_handle)
        end_handle = g.get_previous_step(g.path_end(path_handle))
        if g.get_id(g.get_handle_of_step(start_handle)) in split_node_map_dict and g.get_is_reverse(g.get_handle_of_step(start_handle)):
            patch_handles = bdsg.std.vector_handlegraph_handle_t([g.get_handle(node_id, True) for node_id in split_node_map_dict[g.get_id(g.get_handle_of_step(start_handle))][::-1]])
            g.rewrite_segment(start_handle, g.get_next_step(start_handle), patch_handles)
            patch_i += 1
        if g.get_id(g.get_handle_of_step(end_handle)) in split_node_map_dict and not g.get_is_reverse(g.get_handle_of_step(end_handle)):
            patch_handles = bdsg.std.vector_handlegraph_handle_t([g.get_handle(node_id, False) for node_id in split_node_map_dict[g.get_id(g.get_handle_of_step(end_handle))]])
            g.rewrite_segment(end_handle, g.path_end(path_handle), patch_handles)
            patch_i += 1
    print("{} path ends are patched.".format(patch_i))
    
    
    for ref_path in ref_path_list:
        path_contig = ref_path.split("#")[2]
        if path_contig not in list(vcf_file.header.contigs):
            continue        
        path_range = ref_path.split("#")[3]
        path_handle = g.get_path_handle(ref_path)
        ####Select the reference path region
        if "[" in path_range and "]" in path_range and "-" in path_range:
            path_start = int(path_range[path_range.index("[") + 1 : path_range.index("]")].split("-")[0])
        else:
            path_start = 0
        path_end = path_start
        
        path_step = g.path_begin(path_handle)
        while(g.has_next_step(path_step)):
            path_end += g.get_length(g.get_handle_of_step(path_step))
            path_step = g.get_next_step(path_step)
        path_end += g.get_length(g.get_handle_of_step(path_step))
        

        ####Cycle node
        cycle_node_dict = {}
        path_step = g.path_begin(path_handle)
        while(g.has_next_step(path_step)):
            node_id = g.get_id(g.get_handle_of_step(path_step))
            cycle_node_dict[node_id] = cycle_node_dict.get(node_id, 0) + 1
            path_step = g.get_next_step(path_step)
        node_id = g.get_id(g.get_handle_of_step(path_step))
        cycle_node_dict[node_id] = cycle_node_dict.get(node_id, 0) + 1
        for node, count in cycle_node_dict.items():
            if count >= 2:
                cycle_node_dict[node] = True
            else:
                cycle_node_dict[node] = False
        

        ####Extract the candidate handle paths
        variant_handles_dict = {}
        flanking_handles_dict = {}
        chain_path_list = []
        chain_pos_list = []
        chain_path = []
        node_start = path_start
        node_end = path_start
        path_step = g.path_begin(path_handle)
        variant_i = 0
        while(g.has_next_step(path_step)):
            node_handle = g.get_handle_of_step(path_step)
            node_id = g.get_id(node_handle)
            node_start = node_end
            node_end = node_start + g.get_length(node_handle)
            
            if variant_i < len(variant_list):
                variant_start = variant_list[variant_i][0]
                variant_end = variant_list[variant_i][1]
            else:
                variant_start = path_end
                variant_end = path_end

            if variant_start <= node_start and variant_end >= node_end:
                for variant_record in vcf_file.fetch(path_contig, variant_start, variant_start + 1):
                    if node_start == variant_start:
                        pre_node_start = node_start
                        if variant_record.alts != None:
                            chain_path.append("V" + str(pre_node_start))
                            variant_handles_dict[pre_node_start] = []
                            variant_handles_dict[pre_node_start].append([node_handle])
                            for alt in variant_record.alts:
                                allele_handle = g.create_handle(alt.upper())
                                variant_handles_dict[pre_node_start].append([allele_handle])
                    else:
                        variant_handles_dict[pre_node_start][0].append(node_handle)
                    if node_end == variant_end:
                        variant_i += 1

            else:
                ####start
                if chain_path == []:
                    if not cycle_node_dict[node_id]:
                        pre_node_start = node_start
                        chain_path.append("F" + str(pre_node_start))
                        flanking_handles_dict[pre_node_start] = [node_handle]
                ####variant transition
                elif chain_path[-1][0] == "V":
                    pre_node_start = node_start
                    chain_path.append("F" + str(pre_node_start))
                    flanking_handles_dict[pre_node_start] = [node_handle]
                ####extentsion
                elif cycle_node_dict[node_id]:
                    flanking_handles_dict[pre_node_start].append(node_handle)
                ####end
                else:
                    add_path = False
                    for chain_item in chain_path:
                        if chain_item[0] == "V":
                            add_path = True
                            break
                    if add_path:
                        chain_pos_list.append([int(chain_path[0][1:]), node_end])
                        chain_path_list.append(chain_path)
                        for i, chain_item in enumerate(chain_path):
                            if chain_item[0] == "V":
                                if i != 0:
                                    if chain_path[i - 1][0] == "V":
                                        for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                            for previous_handles in variant_handles_dict[int(chain_path[i - 1][1:])]:
                                                g.create_edge(previous_handles[-1], handles[0])
                                    else:
                                        for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                            previous_handles = flanking_handles_dict[int(chain_path[i - 1][1:])]
                                            g.create_edge(previous_handles[-1], handles[0])
                                if i != len(chain_path) - 1:
                                    if chain_path[i + 1][0] == "V":
                                        for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                            for next_handles in variant_handles_dict[int(chain_path[i + 1][1:])]:
                                                g.create_edge(handles[-1], next_handles[0])
                                    else:
                                        for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                            next_handles = flanking_handles_dict[int(chain_path[i + 1][1:])]
                                            g.create_edge(handles[-1], next_handles[0])
                    ####restart
                    chain_path = []
                    pre_node_start = node_start
                    chain_path.append("F" + str(pre_node_start))
                    flanking_handles_dict[pre_node_start] = [node_handle]
            path_step = g.get_next_step(path_step)
        
        node_handle = g.get_handle_of_step(path_step)
        node_id = g.get_id(node_handle)
        path_step = g.get_next_step(path_step)
        node_start = node_end
        node_end = node_start + g.get_length(node_handle)
        if chain_path != []:
            if chain_path[-1][0] == "V":
                pre_node_start = node_start
                chain_path.append("F" + str(pre_node_start))
                flanking_handles_dict[pre_node_start] = [node_handle]
            else:
                flanking_handles_dict[pre_node_start].append(node_handle)
        
        add_path = False
        for chain_item in chain_path:
            if chain_item[0] == "V":
                add_path = True
                break
        if add_path:
            chain_pos_list.append([int(chain_path[0][1:]), node_end])
            chain_path_list.append(chain_path)
            for i, chain_item in enumerate(chain_path):
                if chain_item[0] == "V":
                    if i != 0:
                        if chain_path[i - 1][0] == "V":
                            for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                for previous_handles in variant_handles_dict[int(chain_path[i - 1][1:])]:
                                    g.create_edge(previous_handles[-1], handles[0])
                        else:
                            for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                previous_handles = flanking_handles_dict[int(chain_path[i - 1][1:])]
                                g.create_edge(previous_handles[-1], handles[0])
                    if i != len(chain_path) - 1:
                        if chain_path[i + 1][0] == "V":
                            for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                for next_handles in variant_handles_dict[int(chain_path[i + 1][1:])]:
                                    g.create_edge(handles[-1], next_handles[0])
                        else:
                            for handles in variant_handles_dict[int(chain_path[i][1:])]:
                                next_handles = flanking_handles_dict[int(chain_path[i + 1][1:])]
                                g.create_edge(handles[-1], next_handles[0])

        print("{} haplotype chains are found.".format(len(chain_path_list)))
        
        ####Add the handle paths
        for chain_path in chain_path_list:
            sample_path_dict = {}
            for sample in sample_list:
                sample_path_dict[sample] = []
                sample_path_dict[sample].append(g.create_path_handle("#".join([sample + args.sample_suffix, "0", args.contig_prefix + str(chain_i), "0"])))
                sample_path_dict[sample].append(g.create_path_handle("#".join([sample + args.sample_suffix, "1", args.contig_prefix + str(chain_i), "0"])))
            for item in chain_path:
                if item[0] == "V":
                    for variant_record in vcf_file.fetch(path_contig, int(item[1:]), int(item[1:]) + 1):
                        for sample in sample_list:
                            gt = (0,0) if variant_record.samples[sample]["GT"] == (None, None) else variant_record.samples[sample]["GT"]
                            for handles in variant_handles_dict[int(item[1:])][gt[0]]:
                                g.append_step(sample_path_dict[sample][0], handles)
                            for handles in variant_handles_dict[int(item[1:])][gt[1]]:
                                g.append_step(sample_path_dict[sample][1], handles)
                else:
                    for sample in sample_list:
                        for handles in flanking_handles_dict[int(item[1:])]:
                            g.append_step(sample_path_dict[sample][0], handles)
                            g.append_step(sample_path_dict[sample][1], handles)
            chain_i += 1
    g.serialize(args.output_graph)


if __name__ == "__main__":
    main()
