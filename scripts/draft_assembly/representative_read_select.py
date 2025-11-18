#!/usr/bin/env python3

# Representative read select script
#
# Select the representative read and calculate the reference coverage proportion of this read.
# 
# Date: 2022/12/10

import sys

def FirstElem(elem):
    return elem[0]

def openfile(filename):
    if filename == "-":
        return sys.stdin
    else:
        return open(filename, "r")

def max_cover_read(bed_file, interval):
    interval_start = int(interval.split(":")[1].split("-")[0])
    interval_end = int(interval.split(":")[1].split("-")[1])

    reads_dict={}
    for line in openfile(bed_file):
        read = line.split()[3]
        start = int(line.split()[1])
        end = int(line.split()[2])
        
        if (start < interval_start):
            start = interval_start
        elif (start > interval_end):
            continue
        if (end > interval_end):
            end = interval_end
        elif (end < interval_start):
            continue

        pos1 = [start, "L"]
        pos2 = [end, "R"]
        if read not in reads_dict.keys():
            reads_dict[read] = [pos1, pos2]
        else:
            reads_dict[read].append(pos1)
            reads_dict[read].append(pos2)

    max_proportion = 0
    max_read = "Nan"
    for read, pos_list in reads_dict.items():
        read_len = 0
        if len(pos_list)==2:
            read_len = pos_list[1][0] - pos_list[0][0] + 1
        else:
            last_pos = 0
            tag = 0
            pos_list.sort(key = FirstElem)
            for pos in pos_list:
                if tag == 0:
                    if last_pos == 0:
                        last_start = pos[0]
                    elif last_pos != pos[0]:
                        last_end = last_pos
                        read_len += last_end - last_start + 1
                        last_start = pos[0]

                if pos[1] == "L":
                    tag+=1
                else:
                    tag-=1
                last_pos = pos[0]
            last_end = last_pos
            read_len += last_end - last_start + 1
        
        proportion = read_len/(interval_end - interval_start +1)
        if(proportion > max_proportion):
            max_proportion = proportion
            max_read = read

    return [max_read, max_proportion]

