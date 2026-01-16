#!/usr/bin/env python3
import sys
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

input_file = sys.argv[1]
final_sample_list = sys.argv[2].split(",")
assembly_sample = sys.argv[3]
sex = sys.argv[4]
if assembly_sample in final_sample_list:
    assembly_index = final_sample_list.index(assembly_sample)
else:
    assembly_index = "NA"

sample_hap_count_list = [0] * len(final_sample_list)
sample_index = []
for line in openfile(input_file):
    if line.startswith("#"):
        sample_list = line.strip().split()[9:]
        for sample in sample_list:
            if sample in final_sample_list:
                sample_index.append(final_sample_list.index(sample))
    if not line.startswith("#"):
        sample_gt_list = line.strip().split()[9:]
        for i, sample_gt in enumerate(sample_gt_list):
            sample_hap_count_list[sample_index[i]] = len(sample_gt.replace("/","|").split("|"))
        break


for line in openfile(input_file):
    if line.startswith("##"):
        print(line.strip())
    elif line.startswith("#"):
        sample_add_list = []
        for i, sample in enumerate(final_sample_list):
            sample_hap_count = sample_hap_count_list[i]
            if sample_hap_count <= 2:
                sample_add_list.append(sample)
            else:
                for j in range(int(sample_hap_count/2)):
                    sample_add_list.append(sample + "-" + str(j + 1))
        print('\t'.join(line.strip().split()[0:9] + sample_add_list))
    else:
        chr = line.strip().split()[0]
        if sex == "female" and "chrY" in chr:
            continue
        sample_gt_list = line.strip().split()[9:]
        sample_gt_add_list = []
        final_sample_gt_list = [".|."] * len(final_sample_list)
        for i, sample_gt in enumerate(sample_gt_list):
            final_sample_gt_list[sample_index[i]] = sample_gt
        for i, sample_gt in enumerate(final_sample_gt_list):
            hap_gt_list = sample_gt.replace("/","|").split("|")
            if len(hap_gt_list) == 1:
                hap_gt_list.append(hap_gt_list[0])
            if sex == "male" and ("chrX" in chr or "chrY" in chr) and i == assembly_index and len(hap_gt_list) == 2:
                if hap_gt_list[0] == "." and hap_gt_list[1] != ".":
                    hap_gt_list[0] = hap_gt_list[1]
                elif hap_gt_list[0] != "." and hap_gt_list[1] == ".":
                    hap_gt_list[1] = hap_gt_list[0]
            for j in range(int(len(hap_gt_list)/2)):
                sample_gt_add_list.append(hap_gt_list[j * 2] + "|" + hap_gt_list[j * 2 + 1])
        print('\t'.join(line.strip().split()[0:9] + sample_gt_add_list))