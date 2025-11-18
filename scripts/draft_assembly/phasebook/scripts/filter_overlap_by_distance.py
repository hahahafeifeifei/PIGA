import numpy as np
import sys

max_distance = 100000
input_file = sys.argv[1]

bed = open(input_file)
bed_dict = {}
for line in bed:
    reads = line.split()[3]
    pos = line.split()[0:3]
    bed_dict[reads] = pos
bed.close()

alignment = sys.stdin
for line in alignment:
    reads1 = line.split()[0]
    reads1_pos = bed_dict[reads1]
    reads2 = line.split()[5]
    reads2_pos = bed_dict[reads2]
    if reads1_pos[0] == reads2_pos[0]:
        distance = np.max([np.max([int(reads1_pos[1]), int(reads2_pos[1])]) - np.min([int(reads1_pos[2]), int(reads2_pos[2])]), 0])
        if distance <= max_distance:
            print(line.strip())
