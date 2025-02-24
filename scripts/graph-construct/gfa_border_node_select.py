import sys
import gzip
import json
import odgi

#Determine and split the boundary nodes of subgraphs in for each about 20 Mb region
#python3 gfa_border_node_select.py input.gfa input.snarls.txt input.og output.node_split.gfa output.node_split.edge_clip.gfa output.bed

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def dist(node_start, node_end, chunk):
    return max([0, node_start - chunk, chunk - node_end])

chunk_size = 20000000
chunk_pos_deviation = 500000

gfa_file = sys.argv[1]
snarls_file = sys.argv[2]
og_file = sys.argv[3]
out_gfa_file = sys.argv[4]
out_clip_edge_gfa_file = sys.argv[5]
out_chunk_file = sys.argv[6]

gfa_dict = {}
paths = {}

# Read the GFA file and store node lengths and paths
for line in openfile(gfa_file):
    if line[0] == "S":
        node_id = line.strip().split()[1]
        node_len = len(line.strip().split()[2])
        gfa_dict[node_id] = node_len
    if line[0] == "P":
        path_name = line.split()[1]
        path_nodes = line.split()[2].replace("+","").replace("-","").split(",")
        paths[path_name] = path_nodes

# Function to process each path
odgi_graph = odgi.graph()
odgi_graph.load(og_file)
all_final_chunks = []
for path_name, ref_path in paths.items():
    ref_node_info = []
    start = 1
    for node in ref_path:
        node_len = gfa_dict[node]
        end = start + node_len - 1
        ref_node_info.append([node, start, end, node_len, True])
        start = end + 1
    ref_len = end

    initial_chunk = []
    # Create initial chunk
    for i in range(int(ref_len/chunk_size)):
        initial_chunk.append([i * chunk_size + 1, (i + 1) * chunk_size])
    initial_chunk.append([(i + 1) * chunk_size + 1, ref_len])

    # Filter in bubble node
    for line in openfile(snarls_file):
        snarls_json = json.loads(line)
        source_node = snarls_json['start']['node_id']
        sink_node = snarls_json['end']['node_id']
        if source_node in ref_path and sink_node in ref_path:
            node_pos = sorted([ref_path.index(source_node), ref_path.index(sink_node)])
            for i in range(node_pos[0] + 1, node_pos[1]):
                ref_node_info[i][4] = False

    # Select the unbubbled nodes and get the number of chunk breakpoints spanned by these nodes
    unbubble_node_info = []
    for node in ref_node_info:
        if node[4]:
            node_span = 0
            for chunk in initial_chunk:
                if dist(node[1], node[2], chunk[0]) == 0:
                    node_span += 1
            unbubble_node_info.append(node[:4] + [node_span])

    # Determine the breakpoint of final chunk based on the distance between selected nodes and initial chunk
    node_min_len = 50000
    final_chunk = []
    for n in range(len(initial_chunk)):
        final_chunk.append([[]] * 2)
    while node_min_len > 0:
        for n, chunk in enumerate(initial_chunk):
            min_dist = [5000000000] * 2
            min_dist_node = [[]] * 2
            for i in [0, 1]:
                for node in unbubble_node_info:
                    if node[3] >= node_min_len:
                        if dist(node[1], node[2], chunk[i]) < min_dist[i]:
                            min_dist[i] = dist(node[1], node[2], chunk[i])
                            min_dist_node[i] = node
                if min_dist[i] <= chunk_pos_deviation and final_chunk[n][i] == []:
                    if min_dist_node[i][4] == 0:
                        node_pos = int((min_dist_node[i][3] + 1)/2) + 1 - i
                        chunk_pos = min_dist_node[i][1] + node_pos - 1
                    elif min_dist_node[i][4] >= 1:
                        if dist(min_dist_node[i][1], min_dist_node[i][2], chunk[i]) == 0:
                            chunk_pos = chunk[i]
                            node_pos = int(chunk[i] - min_dist_node[i][1] + 1)
                        else:
                            continue
                    node_id = min_dist_node[i][0]
                    final_chunk[n][i] = [chunk_pos, node_id, node_pos]
        node_min_len -= 5000
    final_chunk[0][0] = [1, ref_node_info[0][0], 1]
    final_chunk[-1][1] = [ref_len, ref_node_info[-1][0], ref_node_info[-1][3]]

    # Split the chosen node using ODGI API
    node_clip_dist = {}
    for chunk in final_chunk[:-1]:
        node = int(chunk[1][1])
        pos = chunk[1][2]
        if node not in node_clip_dist.keys():
            node_clip_dist[node] = [pos]
        else:
            node_clip_dist[node].append(pos)

    node_map_dist = {}
    for node, pos_list in node_clip_dist.items():
        new_node_handle = []
        node_handle = odgi_graph.get_handle(node, False)
        for pos in pos_list[::-1]:
            split_node_handle = odgi_graph.divide_handle(node_handle, pos)
            node_handle = split_node_handle[0]
            new_node_handle.insert(0, split_node_handle[1])
        new_node_handle.insert(0, node_handle)
        for handle in new_node_handle:
            new_node = odgi_graph.get_id(handle)
            if node not in node_map_dist.keys():
                node_map_dist[node] = [new_node]
            else:
                node_map_dist[node].append(new_node)

    # Add the new node information to the final chunk
    for n in range(len(final_chunk)):
        for i in [0, 1]:
            node = int(final_chunk[n][i][1])
            if node in node_map_dist.keys():
                final_chunk[n][i].append(node_map_dist[node][0])
                if i == 1:
                    del(node_map_dist[node][0])
            else:
                final_chunk[n][i].append(node)
        final_chunk[n].append(path_name)
        all_final_chunks.append(final_chunk[n])


# Write the new GFA file
fo = open(out_gfa_file, "w")
sys.stdout = fo
odgi_graph.to_gfa()
fo.close()

# Write the edge-clipped GFA file
for n in range(len(all_final_chunks) - 1):
    source_node = int(all_final_chunks[n][1][3])
    sink_node = int(all_final_chunks[n + 1][0][3])
    edge = (odgi_graph.get_handle(source_node, False), odgi_graph.get_handle(sink_node, False))
    odgi_graph.destroy_edge(edge)

fo = open(out_clip_edge_gfa_file, "w")
sys.stdout = fo
odgi_graph.to_gfa()
fo.close()

# Write the chunk information for each path
fo = open(out_chunk_file, "w")
for chunk in all_final_chunks:
    fo.write("\t".join(str(i) for i in [chunk[2]] + chunk[0] + chunk[1]) + '\n')
fo.close()