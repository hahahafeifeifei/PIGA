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

def norm_gt_list(gt):
    if gt == ".":
        gt_list = [".", "."]
    else:
        gt1 = gt.replace("|", "/").split("/")[0]
        gt2 = gt.replace("|", "/").split("/")[1]
        if gt1 == ".":
            gt_list = [gt1, gt2]
        elif gt2 == ".":
            gt_list = [gt2, gt1]
        elif int(gt1) <= int(gt2):
            gt_list = [gt1, gt2]
        else:
            gt_list = [gt2, gt1]
    return gt_list

input_vcf_file = sys.argv[1]
output_vcf_file = sys.argv[2]
fo = open(output_vcf_file, "w")
for line in openfile(input_vcf_file):
    if line.startswith("##"):
        fo.write(line)
    elif line.startswith("#"):
        fo.write('\t'.join(line.strip().split()[0:9] + ["sample"]) + '\n')
    else:
        line_list = line.strip().split()
        chr = line_list[0]
        pos = int(line_list[1])
        id = line_list[2]
        ref = line_list[3]
        alt_list = line_list[4].split(",")
        gt_index = line_list[8].split(":").index("GT")
        genotype_gt = line_list[9].split(":")[gt_index]
        genotype_gt_list = norm_gt_list(genotype_gt)
        phase_gt = line_list[10].split(":")[gt_index]
        phase_gt_list = norm_gt_list(phase_gt)
        if genotype_gt_list == phase_gt_list:
            refine_gt = phase_gt
        else:
            refine_gt = genotype_gt
        fo.write('\t'.join(line_list[0:8] + ["GT"] + [refine_gt]) + '\n')
fo.close()
