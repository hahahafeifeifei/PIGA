#!/usr/bin/env python3

# Assembly interval detect script
#
# Detect the assembly intervals with minimum coverage and minimum interval size.
# 
# Date: 2022/12/9

import argparse
import sys
import pysam

def parse_args():
    parser = argparse.ArgumentParser(description='Assembly interval detect script', usage='samtools view -J -Q <min_mq> <haplotagged bam> | assembly_interval_detect.py --haplotagged-bam <haplotagged bam> -o <raw assembly interval list>')
    parser.add_argument('--haplotagged-bam', metavar='<file>', type=str, default='None', help='Input haplotagged bam', required=True)
    parser.add_argument('-o', '--assembly-interval', metavar='<file>', type=str, default='None', help='Output assembly interval list', required=True)
    parser.add_argument('--min-assembly-size', metavar='<int>', type=int, default=10000, help='Minimum interval size of assembly region [10000]')
    parser.add_argument('--min-assembly-cov', metavar='<int>', type=int, default=3, help='Minimum coverage of assembly region [3]')
    parser.add_argument('--min-assembly-mq', metavar='<int>', type=int, default=20, help='Minimum mapping quality of reads [20]')
    parser.add_argument('--min-connect-align-len', metavar='<int>', type=int, default=2000, help='Minimum alignment length used to connect assembly interval [2000]')
    return parser.parse_args()

def raw_assembly_interval_detect(min_assembly_size, min_assembly_cov):
    raw_assembly_interval_list = []
    record_chr = "0"
    record_pos = 0
    pre_chr = "0"
    pre_pos = 0
    for line in sys.stdin:
        depth = int(line.strip().split("\t")[2])
        if depth >= min_assembly_cov:
            chr = line.strip().split("\t")[0]
            pos = int(line.strip().split("\t")[1])
            if chr != pre_chr or pos != pre_pos+1:
                if record_chr!= "0" and pre_pos - record_pos + 1 >= min_assembly_size:
                    raw_assembly_interval_list.append([record_chr, record_pos, pre_pos])
                record_chr = chr
                record_pos = pos
            pre_chr = chr
            pre_pos = pos
    if record_chr!= "0" and pre_pos - record_pos + 1 >= min_assembly_size:
        raw_assembly_interval_list.append([record_chr, record_pos, pre_pos])
    return raw_assembly_interval_list


def interval_overlap(start1, end1, start2, end2):
    max_start = start1 if start1 >= start2 else start2
    min_end = end1 if end1 <= end2 else end2
    overlap_len = min_end - max_start + 1
    return overlap_len if overlap_len >=0 else 0

def clip_reads_link(haplotagged_bam, raw_assembly_interval_list, assembly_interval, min_assembly_cov, min_connect_align_len, min_assembly_mq):
    #load the reads with supplementary alignment
    reads_dict = {}
    for reads in pysam.AlignmentFile(haplotagged_bam,"rb"):
        if not reads.is_secondary and reads.has_tag("SA"):
            if reads.query_alignment_length >= min_connect_align_len and reads.mapping_quality >= min_assembly_mq:
                reads_name = reads.query_name
                align_chr = reads.reference_name
                align_start = reads.reference_start + 1
                align_end = reads.reference_end
                if reads_name not in reads_dict.keys():
                    reads_dict[reads_name] = [[align_chr, align_start, align_end]]
                else:
                    reads_dict[reads_name].append([align_chr, align_start, align_end])
    
    #filter reads with only 1 alignment
    reads_filter_dict = {}
    for reads_name,reads in reads_dict.items():
        if len(reads)!=1:
            reads_filter_dict[reads_name] = reads

    #load raw assembly interval
    interval_chr_dict = {}
    for interval in raw_assembly_interval_list:
        interval_chr = interval[0]
        interval_pos = [interval[1], interval[2], 0]
        if interval_chr not in interval_chr_dict.keys():
            interval_chr_dict[interval_chr] = [interval_pos]
        else:
            interval_chr_dict[interval_chr].append(interval_pos)
    
    #find the interval linked by clip reads
    interval_pair=[]
    for reads_name,reads in reads_filter_dict.items():
        for i in range(len(reads)-1):
            for j in range(i+1,len(reads)):
                if reads[i][0] == reads[j][0]:
                    align_chr = reads[i][0]
                    if align_chr not in interval_chr_dict.keys():
                        continue
                    
                    #find intervals overlapped with alignment1
                    align_start1 = reads[i][1]
                    align_end1 = reads[i][2]
                    is_overlap = False
                    overlap_interval1 = []
                    for interval_n in range(len(interval_chr_dict[align_chr])):
                        interval_start = interval_chr_dict[align_chr][interval_n][0]
                        interval_end = interval_chr_dict[align_chr][interval_n][1]
                        overlap_len = interval_overlap(align_start1, align_end1, interval_start, interval_end)
                        if overlap_len != 0:
                            is_overlap = True
                            overlap_interval1.append([interval_start, interval_end, overlap_len, interval_n])
                        elif overlap_len == 0:
                            if is_overlap:
                                break
                    
                    #find intervals overlapped with alignment2
                    align_start2 = reads[j][1]
                    align_end2 = reads[j][2]
                    is_overlap = False
                    overlap_interval2 = []
                    for interval_n in range(len(interval_chr_dict[align_chr])):
                        interval_start = interval_chr_dict[align_chr][interval_n][0]
                        interval_end = interval_chr_dict[align_chr][interval_n][1]
                        overlap_len = interval_overlap(align_start2, align_end2, interval_start, interval_end)
                        if overlap_len != 0:
                            is_overlap = True
                            overlap_interval2.append([interval_start, interval_end, overlap_len, interval_n])
                        elif overlap_len == 0:
                            if is_overlap:
                                break
                    
                    #find the largest overlap interval pair with alignments
                    overlap_len_sum = 0
                    interval_pair_info = []
                    for interval1 in overlap_interval1:
                        for interval2 in overlap_interval2:
                            if abs(interval1[3] - interval2[3]) == 1 and (interval1[2] + interval2[2]) > overlap_len_sum:
                                overlap_len_sum = interval1[2] + interval2[2]
                                if interval1[0] <= interval2[0]:
                                    interval_pair_info = [reads_name, align_chr, interval1[0], interval1[1], interval1[3], interval2[0], interval2[1], interval2[3]]
                                else:
                                    interval_pair_info = [reads_name, align_chr, interval2[0], interval2[1], interval2[3], interval1[0], interval1[1], interval1[3]]

                    #update the number of supporting reads
                    if interval_pair_info != [] and interval_pair_info not in interval_pair:
                        interval_pair.append(interval_pair_info) 
                        interval_chr_dict[align_chr][interval_pair_info[4]][2] += 1

    fo = open(assembly_interval, "w")
    for interval_chr, interval_chr_info in interval_chr_dict.items():
        merge = False
        for interval_info in interval_chr_info:
            if not merge:
                fo.write(interval_chr + '\t' + str(interval_info[0]) + "\t")
            if interval_info[2] >= min_assembly_cov:
                merge = True
            else:
                merge = False
            if not merge:
                fo.write(str(interval_info[1]) + '\n')
    fo.close()

def assembly_interval_split(assembly_interval_list, max_sub_interval_size, sub_interval_overlap_size):
    assembly_split_interval_dict = {}
    for line in open(assembly_interval_list, "r"):
        chr = line.split()[0]
        start = int(line.split()[1])
        end = int(line.split()[2])
        interval_size = end - start + 1
        assembly_split_interval_dict[chr + ':' + str(start) + "-" + str(end)] = []
        if interval_size > max_sub_interval_size:
            sub_interval_size = interval_size
            sub_start = start
            while sub_interval_size > max_sub_interval_size:
                sub_end = sub_start + max_sub_interval_size - 1
                assembly_split_interval_dict[chr + ':' + str(start) + "-" + str(end)].append(chr + ':' + str(sub_start) + "-" + str(sub_end))
                sub_start = sub_start + sub_interval_overlap_size
                sub_interval_size = end - sub_start + 1
            assembly_split_interval_dict[chr + ':' + str(start) + "-" + str(end)].append(chr + ':' + str(sub_start) + "-" + str(end))
        else:
            assembly_split_interval_dict[chr + ':' + str(start) + "-" + str(end)].append(chr + ':' + str(start) + "-" + str(end))
    return assembly_split_interval_dict

def main():
    args = parse_args()
    raw_assembly_interval_list = raw_assembly_interval_detect(args.min_assembly_size, args.min_assembly_cov)
    clip_reads_link(args.haplotagged_bam, raw_assembly_interval_list, args.assembly_interval, args.min_assembly_cov, args.min_connect_align_len, args.min_assembly_mq)

if __name__ == "__main__":
    main()

