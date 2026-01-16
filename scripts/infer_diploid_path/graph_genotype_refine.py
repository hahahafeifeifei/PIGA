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
pangenie_sample = sys.argv[3]
snp_sample = sys.argv[4]
assembly_sample = sys.argv[5]

fo = open(output_vcf_file, "w")
for line in openfile(input_vcf_file):
    if line.startswith("##"):
        fo.write(line)
    elif line.startswith("#"):
        fo.write('\t'.join(line.strip().split()[0:9] + ["sample"]) + '\n')
        pangenie_index = line.strip().split().index(pangenie_sample)
        snp_index = line.strip().split().index(snp_sample)
        assembly_index = line.strip().split().index(assembly_sample)
    else:
        line_list = line.strip().split()
        chr = line_list[0]
        pos = int(line_list[1])
        id = line_list[2]
        ref = line_list[3]
        alt_list = line_list[4].split(",")
        if "GT" not in line_list[8].split(":") or "GQ" not in line_list[8].split(":"):
            continue
        gt_index = line_list[8].split(":").index("GT")
        gq_index = line_list[8].split(":").index("GQ")
        gq = 0 if line_list[pangenie_index].split(":")[gq_index] == "." else int(line_list[pangenie_index].split(":")[gq_index])
        pangenie_gt_list = norm_gt_list(line_list[pangenie_index].split(":")[gt_index])
        snp_gt_list = norm_gt_list(line_list[snp_index].split(":")[gt_index])
        assembly_gt_list = norm_gt_list(line_list[assembly_index].split(":")[gt_index])
        if gq <= 50:
            if len(alt_list) == 1 and "." not in snp_gt_list:
                refine_gt_list = snp_gt_list
            elif "." not in assembly_gt_list:
                refine_gt_list = assembly_gt_list
            else:
                refine_gt_list = pangenie_gt_list
        else:
            refine_gt_list = pangenie_gt_list
        refine_gt = "/".join(refine_gt_list)
        fo.write('\t'.join(line_list[0:8] + ["GT"] + [refine_gt]) + '\n')
fo.close()