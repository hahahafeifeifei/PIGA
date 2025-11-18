#!/usr/bin/env python3

# Pacbio zmw haplotag script
#
# Partition pacbio zmws into two haplotypes based on the whatshap haplotag results of hifi reads and subreads.
# We choose the majority voting haplotype as subreads haplotype. And we use the whatshap-partitioned haplotype as hifi reads haplotype.
# Finally, we generate the union haplotype set of subreads and hifi reads, trusting hifi reads haplotype more.
# 
# Date: 2022/12/6

import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Pacbio zmw haplotag script', usage='pacbio_zmw_haplotag.py --subread-haplotag <subreads haplotag list> --hifi-haplotag <hifi reads haplotag list> -o <consensus haplotag list>')
    parser.add_argument('--subread-haplotag', metavar='<file>', type=str, default='None', help='Input Pacbio subread haplotag file', required=True )
    parser.add_argument('--hifi-haplotag', metavar='<file>', type=str, default='None', help='Input Pacbio hifi read haplotag file', required=True)
    parser.add_argument('-o', '--consensus-haplotag', metavar='<file>', type=str, default='None', help='Output consensus haplotag file', required=True)
    return parser.parse_args()

def openfile(filename):
    if filename == "-":
        return sys.stdin
    else:
        return open(filename, "r")

def consensus_haplotag(subreads_haplotag_list, hifi_haplotag_list, consensus_haplotag_list):
    zmw_hap_dict = {}
    subreads_zmw_hap_dict = {}
    subreads_zmw_reads_dict = {}
    for line in openfile(subreads_haplotag_list):
        if(line[0] == "#"):
            continue
        read = line.strip().split()[0]
        zmw = '/'.join(read.split('/')[0:2])
        hap = line.strip().split()[1]

        if zmw not in subreads_zmw_hap_dict.keys():
            subreads_zmw_hap_dict[zmw] = [0, 0, 0]
            subreads_zmw_reads_dict[zmw] = [read]
        else:
            subreads_zmw_reads_dict[zmw].append(read)

        if hap == "H1":
            subreads_zmw_hap_dict[zmw][0] += 1
        elif hap == "H2":
            subreads_zmw_hap_dict[zmw][1] += 1
        else:
            subreads_zmw_hap_dict[zmw][2] += 1

    for zmw, zmw_list in subreads_zmw_hap_dict.items():
        if zmw_list[2]/(zmw_list[0] + zmw_list[1] + zmw_list[2]) > 0.5:
            zmw_hap_dict[zmw] = ["none","none"]
        elif zmw_list[0] > zmw_list[1]:
            zmw_hap_dict[zmw] = ["H1","none"]
        elif zmw_list[0] < zmw_list[1]:
            zmw_hap_dict[zmw] = ["H2","none"]
        elif zmw_list[0] == zmw_list[1]:
            zmw_hap_dict[zmw] = ["none","none"]

    for line in openfile(hifi_haplotag_list):
        if(line[0] == "#"):
            continue
        read = line.strip().split()[0]
        zmw = '/'.join(read.split('/')[0:2])
        hap = line.strip().split()[1] 

        if zmw not in zmw_hap_dict.keys():
            zmw_hap_dict[zmw] = ['none', hap]
        else:
            zmw_hap_dict[zmw][1] = hap

    fo = open(consensus_haplotag_list, 'w')
    for zmw, hap_list in zmw_hap_dict.items():
        if hap_list[0] == hap_list[1]:
            hap = hap_list[0]
        elif hap_list[1] == "none":
            hap = hap_list[0]
        else:
            hap = hap_list[1]
        
        if zmw in subreads_zmw_reads_dict.keys():
            for read in subreads_zmw_reads_dict[zmw]:
                fo.write(read + '\t' + hap + '\t' + zmw + '\n')
    fo.close()
        

def main():
    args = parse_args()
    consensus_haplotag(args.subread_haplotag, args.hifi_haplotag, args.consensus_haplotag)

if __name__ == "__main__":
    main()
