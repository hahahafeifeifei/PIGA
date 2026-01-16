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

def gt_list(gt):
    if gt == ".":
        gt_list = [".", "."]
    else:
        gt_list = gt.replace("|", "/").split("/")
    return gt_list

def norm_gt_list(gt):
    if gt == ".":
        norm_gt_list = [".", "."]
    else:
        gt1 = gt.replace("|", "/").split("/")[0]
        gt2 = gt.replace("|", "/").split("/")[1]
        if gt1 == ".":
            norm_gt_list = [gt1, gt2]
        elif gt2 == ".":
            norm_gt_list = [gt2, gt1]
        elif int(gt1) <= int(gt2):
            norm_gt_list = [gt1, gt2]
        else:
            norm_gt_list = [gt2, gt1]
    return norm_gt_list

input_vcf_file = sys.argv[1]
output_vcf_file = sys.argv[2]
target_sample = sys.argv[3]
phase_rank_sample_list = sys.argv[4].split(",")

phase_index_list = []
pre_chr = ""
fo = open(output_vcf_file, "w")
for line in openfile(input_vcf_file):
    if line.startswith("##"):
        fo.write(line)
    elif line.startswith("#"):
        info_line = line.strip().split()
        for sample in phase_rank_sample_list:
            phase_index_list.append(info_line.index(sample))
        target_index = info_line.index(target_sample)
        fo.write('\t'.join(line.strip().split()[0:9] + ["sample"]) + '\n')
    else:
        line_list = line.strip().split()
        chr = line_list[0]
        pos = int(line_list[1])
        id = line_list[2]
        ref = line_list[3]
        alt_list = line_list[4].split(",")
        gt_index = line_list[8].split(":").index("GT")
        is_unphase = False

        if chr != pre_chr:
            pre_chr = chr
            phase_filp_list = [None] * len(phase_index_list)

        target_gt = line_list[target_index].split(":")[gt_index]
        target_gt_list = gt_list(target_gt)
        if target_gt_list[0] != target_gt_list[1]:
            if phase_filp_list != [None] * len(phase_index_list):
                ####target phased genotype determination
                cis_list = []
                trans_list = []
                for i, phase_index in enumerate(phase_index_list):
                    phase_gt = line_list[phase_index].split(":")[gt_index]
                    phase_gt_list = gt_list(phase_gt)
                    if phase_filp_list[i] == None:
                        continue
                    elif phase_filp_list[i] == True:
                        phase_gt_list.reverse()
                    #print(target_gt_list)
                    if "|" in phase_gt:
                        if target_gt_list[0] == phase_gt_list[0] and target_gt_list[1] == phase_gt_list[1]:
                            cis_list.append(i)
                        if target_gt_list[1] == phase_gt_list[0] and target_gt_list[0] == phase_gt_list[1]:
                            trans_list.append(i)
                #print(phase_filp_list,end =" ")
                print(cis_list, end=" ")
                print(trans_list)

                if len(cis_list) == len(trans_list):
                    if len(cis_list) == 0:
                        is_unphase = True
                    else:
                        if cis_list[0] > trans_list[0]:
                            target_gt_list.reverse()
                elif len(cis_list) < len(trans_list):
                    target_gt_list.reverse()
            #print(phase_filp_list, end=" ")
            ####phase phased genotype flipping
            if not is_unphase:
                for i, phase_index in enumerate(phase_index_list):
                    phase_gt = line_list[phase_index].split(":")[gt_index]
                    phase_gt_list = gt_list(phase_gt)
                    #print(target_gt_list, end=" ")
                    #print(phase_gt_list, end=" ")
                    if "|" in phase_gt:
                        if target_gt_list[0] == phase_gt_list[0] and target_gt_list[1] == phase_gt_list[1]:
                            phase_filp_list[i] = False
                        if target_gt_list[1] == phase_gt_list[0] and target_gt_list[0] == phase_gt_list[1]:
                            phase_filp_list[i] = True
                    #print(phase_filp_list[i])
            #print(phase_filp_list)
        if is_unphase:
            target_gt = "/".join(target_gt_list)
        else:
            target_gt = "|".join(target_gt_list)
        fo.write('\t'.join(line_list[0:8] + ["GT"] + [target_gt]) + '\n')
fo.close()