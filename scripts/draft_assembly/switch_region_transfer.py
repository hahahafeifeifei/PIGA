#!/usr/bin/env python3
import pysam
import argparse

# Switch region transfer script
#
# Transfer the clip region from reference to assembly.
# 
# Date: 2022/12/16

def parse_args():
    parser = argparse.ArgumentParser(description='Switch region transfer script', usage='switch_region_transfer.py -c <input bed> -b <alignment bam> -o <output bed>')
    parser.add_argument('-c', '--clip-ref-bed', metavar='<file>', type=str, default='None', help='Input bed file of switch region in reference', required=True)
    parser.add_argument('-b', '--contig-bam', metavar='<file>', type=str, default='None', help='Input contig alignment bam file', required=True)
    parser.add_argument('-o', '--clip-assembly-bed', metavar='<file>', type=str, default='None', help='Output bed file of switch region in assembly', required=True)
    return parser.parse_args()

def openfile(filename):
    if filename.endswith('.bam'):
        return pysam.AlignmentFile(filename,"rb")
    else:
        return open(filename, "r")

def transfer(args):
    bed_list = []
    bam_file = openfile(args.contig_bam)
    for region in openfile(args.clip_ref_bed):
        interval_chr = region.strip().split()[0]
        interval_start = int(region.strip().split()[1])
        interval_end = int(region.strip().split()[2])
        start_alignment_list = []
        end_alignment_list = []
        all_alignment_list = []

        for base in bam_file.pileup(contig=interval_chr, start=interval_start, end=interval_end, truncate=True):
            if base.reference_pos == interval_start:
                for start_read in base.pileups:
                    start_alignment_list.append(start_read)

            elif base.reference_pos == interval_end - 1:
                for end_read in base.pileups:
                    end_alignment_list.append(end_read)

            for read in base.pileups:
                if read.alignment not in all_alignment_list:
                    all_alignment_list.append(read.alignment)
            
        for read in all_alignment_list:
            read_name = read.query_name
            read_start = read.query_alignment_start
            read_end = read.query_alignment_end
            
            if read.reference_start < interval_start:
                for start_read in start_alignment_list:
                    if read == start_read.alignment:
                        if start_read.is_del:
                            read_start = start_read.query_position_or_next
                        else:
                            read_start = start_read.query_position
                        break
             
            if read.reference_end > interval_end:
                for end_read in end_alignment_list:
                    if read == end_read.alignment:
                        if end_read.is_del:
                            read_end = end_read.query_position_or_next
                        else:
                            read_end = end_read.query_position
                        break            
            if str(read_start) != None and str(read_end) != None:
                bed_info = [read_name, str(read_start), str(read_end)]
                bed_list.append(bed_info)

    fo = open(args.clip_assembly_bed, "w")
    for bed_info in bed_list:
        fo.write("\t".join(bed_info)+'\n')
    fo.close()


def main():
    args = parse_args()
    transfer(args)

if __name__ == "__main__":
    main()


