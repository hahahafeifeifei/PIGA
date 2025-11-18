#!/usr/bin/env python3

# Pacbio male extra haplotag script
#
# If the sample is male, after getting the consensus haplotag result, partition extra pacbio zmws on chrX non-PAR region and chrY into two haplotypes based on the aligment results of hifi reads and subreads.
# 
# Date: 2022/12/6

import argparse
import sys
import pysam

def parse_args():
    parser = argparse.ArgumentParser(description='Pacbio male extra haplotag script', usage='pacbio_male_extra_haplotag.py --raw-haplotag <raw haplotag list> --subreads-bam <subread alignment bam> --hifi-bam <hifi read alignment bam> --par-bed <PAR region bed> -o <extra haplotag list>')
    parser.add_argument('--raw-haplotag', metavar='<file>', type=str, default='None', help='Input raw haplotag file', required=True)
    parser.add_argument('--subreads-bam', metavar='<file>', type=str, default='None', help='Input Pacbio subread alignment bam file', required=True)
    parser.add_argument('--hifi-bam', metavar='<file>', type=str, default='None', help='Input Pacbio hifi read alignment bam file', required=True)
    parser.add_argument('--par-bed', metavar='<file>', type=str, default='None', help='Input the PAR region bed', required=True)
    parser.add_argument('-o', '--extra-haplotag', metavar='<file>', type=str, default='None', help='Output extra haplotag file', required=True)
    return parser.parse_args()

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.bam'):
        return pysam.AlignmentFile(filename,"rb")
    elif filename.endswith('.sam'):
        return pysam.AlignmentFile(filename,"r")
    else:
        return open(filename, "r")

def extra_haplotag(raw_haplotag_list, subreads_bam, hifi_bam, par_bed, extra_haplotag_list):
    zmw_reads_hap_chr_dict = {}
    zmw_hap_dict = {}
    for line in openfile(raw_haplotag_list):
        read = line.strip().split()[0]
        hap = line.strip().split()[1]
        zmw = line.strip().split()[2]

        if zmw not in zmw_reads_hap_chr_dict.keys():
            zmw_reads_hap_chr_dict[zmw] = {}
        zmw_reads_hap_chr_dict[zmw][read] = [hap, "A"]  
        if zmw not in zmw_hap_dict.keys():
            zmw_hap_dict[zmw] = hap

    par_info = openfile(par_bed).readlines()
    nonpar_start = int(par_info[0].split("\t")[2]) + 1
    nonpar_end = int(par_info[1].split("\t")[1])
    
    subreads_zmw_chr_dict = {}
    subread_bam_info = openfile(subreads_bam)
    chrx_read = subread_bam_info.fetch("chrX", nonpar_start, nonpar_end)
    for read_info in chrx_read:
        if (not read_info.is_secondary) and (not read_info.is_supplementary):
            read = read_info.query_name
            zmw = '/'.join(read_info.query_name.split('/')[0:2])
            if zmw in zmw_hap_dict.keys() and zmw_hap_dict[zmw] == "none":
                if read in zmw_reads_hap_chr_dict[zmw].keys():
                    zmw_reads_hap_chr_dict[zmw][read][1] = "X"
                    if zmw not in subreads_zmw_chr_dict.keys():
                        subreads_zmw_chr_dict[zmw] = [0, 0, 0]
                    subreads_zmw_chr_dict[zmw][1] += 1
            
    chry_read = subread_bam_info.fetch("chrY")
    for read_info in chry_read:
        if (not read_info.is_secondary) and (not read_info.is_supplementary):
            read = read_info.query_name
            zmw = '/'.join(read_info.query_name.split('/')[0:2])
            if zmw in zmw_hap_dict.keys() and zmw_hap_dict[zmw] == "none":
                if read in zmw_reads_hap_chr_dict[zmw].keys():
                    zmw_reads_hap_chr_dict[zmw][read][1] = "Y"
                    if zmw not in subreads_zmw_chr_dict.keys():
                        subreads_zmw_chr_dict[zmw] = [0, 0, 0]
                    subreads_zmw_chr_dict[zmw][2] += 1

    zmw_chr_dict = {}
    for zmw in subreads_zmw_chr_dict.keys():
        subreads_zmw_chr_dict[zmw][0] = len(zmw_reads_hap_chr_dict[zmw]) - subreads_zmw_chr_dict[zmw][1] - subreads_zmw_chr_dict[zmw][2]
        if subreads_zmw_chr_dict[zmw][1]/(subreads_zmw_chr_dict[zmw][0] + subreads_zmw_chr_dict[zmw][1] + subreads_zmw_chr_dict[zmw][2]) > 0.5:
            zmw_chr_dict[zmw] = ["X", "A"]
        elif subreads_zmw_chr_dict[zmw][2]/(subreads_zmw_chr_dict[zmw][0] + subreads_zmw_chr_dict[zmw][1] + subreads_zmw_chr_dict[zmw][2]) > 0.5:
            zmw_chr_dict[zmw] = ["Y", "A"]

    hifi_bam_info = openfile(hifi_bam)
    chrx_read = hifi_bam_info.fetch("chrX", nonpar_start, nonpar_end)
    for read_info in chrx_read:
        if (not read_info.is_secondary) and (not read_info.is_supplementary):
            read = read_info.query_name
            zmw = '/'.join(read_info.query_name.split('/')[0:2])
            if zmw in zmw_chr_dict.keys():
                zmw_chr_dict[zmw][1] = "X"

    chry_read = hifi_bam_info.fetch("chrY")
    for read_info in chry_read:
        if (not read_info.is_secondary) and (not read_info.is_supplementary):
            read = read_info.query_name
            zmw = '/'.join(read_info.query_name.split('/')[0:2])
            if zmw in zmw_chr_dict.keys():
                zmw_chr_dict[zmw][1] = "Y"

    for zmw, chr_list in zmw_chr_dict.items():
        if chr_list[0] == chr_list[1]:
            chr = chr_list[0]
        elif chr_list[1] == "A":
            chr = chr_list[0]
        else:
            chr = chr_list[1]
        
        for read in zmw_reads_hap_chr_dict[zmw].keys():
            if chr == "X" and zmw_reads_hap_chr_dict[zmw][read][1] == "X":
                zmw_reads_hap_chr_dict[zmw][read][0] = "H2"
            elif chr == "Y" and zmw_reads_hap_chr_dict[zmw][read][1] == "Y":
                zmw_reads_hap_chr_dict[zmw][read][0] = "H1"

    fo = open(extra_haplotag_list, 'w')
    for zmw in zmw_reads_hap_chr_dict.keys():
        for read in zmw_reads_hap_chr_dict[zmw].keys():
            hap = zmw_reads_hap_chr_dict[zmw][read][0]
            fo.write(read + '\t' + hap + '\t' + zmw + '\n')
    fo.close()

def main():
    args = parse_args()
    extra_haplotag(args.raw_haplotag, args.subreads_bam, args.hifi_bam, args.par_bed, args.extra_haplotag)

if __name__ == "__main__":
    main()
