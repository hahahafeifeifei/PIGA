#!/usr/bin/env python3

# Low-coverage phase-assembly pipeline
# A pipeline to assemble haplotype-resolved genome using low-coverage long-reads (and short-reads) data.
#
# Author: Yifei Wang
# Contact: wangyifei05@westlake.edu.cn
# Version: 2.0.5
# Date: 2022/12/18

import argparse
from configparser import RawConfigParser
import os
import sys
from config_logging import logger
from distutils import spawn
import subprocess
from multiprocessing import Pool
import pysam
from assembly_interval_detect import assembly_interval_split
from representative_read_select import max_cover_read

def parse_args():
    parser = argparse.ArgumentParser(description='phase-assembly pipeline version2.', usage='phase_assembly.version2.py [option]')
    parser.add_argument('--subread-bam', metavar='<file>', type=str, default='None', help='Input Pacbio subread alignment bam file', required=True )
    parser.add_argument('--hifi-bam', metavar='<file>', type=str, default='None', help='Input Pacbio hifi read alignment bam file', required=True)
    parser.add_argument('--ngs-bam', metavar='<file>', type=str, default='None', help='Input NGS read alignment bam file')
    parser.add_argument('-v','--vcf', metavar='<file>', type=str, default='None', help='Input vcf file with phased heterzygous SNPs', required=True)
    parser.add_argument('-r', '--reference', metavar='<file>', type=str, default='None', help='Reference genome fasta file', required=True)
    parser.add_argument('-p','--prefix', metavar='<str>', type=str, default='None', help='Prefix of output', required=True)
    parser.add_argument('-o', '--out-dir', metavar='<path>', type=str, default='./assembly-out/', help='Output directory [./assembly-out/]')
    parser.add_argument('-T','--tmp-dir', metavar='<path>', type=str, default='/tmp', help='Temporary directory [./tmp]')
    parser.add_argument('-g','--gender', metavar='<str>', type=str, default='female', help='The gender of the sample [female]')
    parser.add_argument('--par-bed', metavar='<file>', type=str, default='None', help='The PAR region bed', required=True)
    parser.add_argument('--min-assembly-mq', metavar='<int>', type=int, default=20, help='Minimum mapping quality of reads [20]')
    parser.add_argument('--min-assembly-cov', metavar='<int>', type=int, default=3, help='Minimum coverage of assembly region [3]')
    parser.add_argument('--min-assembly-size', metavar='<int>', type=int, default=10000, help='Minimum interval size of assembly region [10000]')
    parser.add_argument('--min-connect-align-len', metavar='<int>', type=int, default=2000, help='Minimum alignment length used to connect assembly interval [2000]')
    parser.add_argument('--max-sub-interval-size', metavar='<int>', type=int, default=90000, help='Maximum sub interval size of assembly region [90000]')
    parser.add_argument('--sub-interval-overlap-size', metavar='<int>', type=int, default=30000, help='The overlap region size between adjacent region [30000]')
    parser.add_argument('-t', '--threads', metavar='<int>', type=int, default=1, help='Number of threads [1]')
    return parser.parse_args()

def check_software():
    SAMTOOLS = ['samtools', 'bgzip', 'tabix']
    WHATSHAP = ['whatshap']
    PBMM = ['pbmm2']
    MINIMAP = ['minimap2']
    MECAT = ['mecat2map', 'mecat2pm4', 'mecat2lcr', 'mecat2splitreads', 'mecat2trimbases']
    SEQKIT = ['seqkit']
    WTDBG = ['wtdbg2', 'wtpoa-cns']
    RACON = ['racon']
    ARROW = ['gcpp']
    CUTADAPT = ['cutadapt']
    SNIFFLES = ['sniffles']
    BCFTOOLS = ['bcftools']
    BEDTOOLS = ['bedtools']
    JASMINE = ['iris']
    GATK = ['gatk']
    FREEBAYES = ['freebayes']

    for software in [SAMTOOLS, WHATSHAP, PBMM, MINIMAP, MECAT, SEQKIT, WTDBG, ARROW, RACON, CUTADAPT, SNIFFLES, BCFTOOLS, BEDTOOLS, JASMINE, GATK, FREEBAYES]:
        for sub_software in software:
            if not spawn.find_executable(sub_software):
                logger.error(sub_software + " is not installed." )
                sys.exit(1)

def check_file(file):
    if not os.path.isfile(file):
        logger.error("Missing file: " + file)
        sys.exit(1)

def run_command(cmd, max_time = 86400):
    logger.info("Running the command: \"%s\"", cmd)
    run_result = subprocess.run(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, timeout = max_time)
    if run_result.returncode == 0:
        logger.info("Successfully ran command: \"%s\"", cmd)
    else: 
        logger.error("Command ran failed: \"%s\"" + "\n" + "\n" + "Error message: " + "\n" + "%s" + "%s" + '\n' , cmd, run_result.stdout.decode('UTF-8').strip(), run_result.stderr.decode('UTF-8').strip() )    
        raise Exception("Subprocess running error")

def run_command_no_info(cmd, max_time = 86400):
    run_result = subprocess.run(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE, timeout = max_time)
    if run_result.returncode != 0:
        logger.error("Command ran failed: \"%s\"" + "\n" + "\n" + "Error message: " + "\n" + "%s" + "%s" + '\n' , cmd, run_result.stdout.decode('UTF-8').strip(), run_result.stderr.decode('UTF-8').strip() )    
        raise Exception("Subprocess running error")

def haplotag_parallel(chrom, ref, haplotag_list, vcf, bam):
    run_command("whatshap haplotag --ignore-read-groups --regions {} --reference {} --output-haplotag-list {} {} {} >/dev/null".format(chrom, ref, haplotag_list, vcf, bam))

def haplotag(args):
    chroms = list(pysam.VariantFile(args.vcf).header.contigs)
    #######Subread haplotag
    subread_haplotag_pool = Pool(args.threads)
    for chrom in chroms:
        chrom_haplotag_list = os.path.join(args.tmp_dir, args.prefix + "." + chrom + ".subread.haplotag.list")
        subread_haplotag_pool.apply_async(haplotag_parallel, args = (chrom, args.reference, chrom_haplotag_list, args.vcf, args.subread_bam))
    subread_haplotag_pool.close()
    subread_haplotag_pool.join()

    subreads_haplotag_list = os.path.join(args.out_dir, args.prefix + ".subread.haplotag.list")
    run_command("cat {} | grep -v \"^#\" > {}".format(os.path.join(args.tmp_dir, args.prefix + ".*.subread.haplotag.list"), subreads_haplotag_list))
    run_command("rm {}".format(os.path.join(args.tmp_dir, args.prefix + ".*.subread.haplotag.list")))
    check_file(subreads_haplotag_list)    

    #######HiFi read haplotag
    hifi_haplotag_pool = Pool(args.threads)
    for chrom in chroms:
        chrom_haplotag_list = os.path.join(args.tmp_dir, args.prefix + "." + chrom + ".hifi.haplotag.list")
        hifi_haplotag_pool.apply_async(haplotag_parallel, args = (chrom, args.reference, chrom_haplotag_list, args.vcf, args.hifi_bam))
    hifi_haplotag_pool.close()
    hifi_haplotag_pool.join()

    hifi_haplotag_list = os.path.join(args.out_dir, args.prefix + ".hifi.haplotag.list")
    run_command("cat {} | grep -v \"^#\" > {}".format(os.path.join(args.tmp_dir, args.prefix + ".*.hifi.haplotag.list"), hifi_haplotag_list))
    run_command("rm {}".format(os.path.join(args.tmp_dir, args.prefix + ".*.hifi.haplotag.list")))
    check_file(hifi_haplotag_list)

    #######Consensus haplotag
    raw_haplotag_list = os.path.join(args.out_dir, args.prefix + ".raw.haplotag.list")
    run_command("{}/pacbio_zmw_haplotag.py --subread-haplotag {} --hifi-haplotag {} -o {}".format(sys.path[0], subreads_haplotag_list, hifi_haplotag_list, raw_haplotag_list))
    check_file(raw_haplotag_list)

    consensus_haplotag_list = os.path.join(args.out_dir, args.prefix + ".consensus.haplotag.list")
    if args.gender == "male":
        run_command("{}/pacbio_male_extra_haplotag.py --raw-haplotag {} --subreads-bam {} --hifi-bam {} --par-bed {} -o {} ".format(sys.path[0], raw_haplotag_list, args.subread_bam, args.hifi_bam, args.par_bed, consensus_haplotag_list))
    else:
         run_command("cp {} {}".format(raw_haplotag_list, consensus_haplotag_list))
    run_command("rm {}".format(raw_haplotag_list))
    check_file(consensus_haplotag_list)

def freebayes_parallel(bam, ref, chrom, vcf):
    run_command("freebayes -b {} -f {} -r {} | bgzip -c > {}".format(bam, ref, chrom, vcf))
    run_command("tabix -f {}".format(vcf))

def interval_assembly_parallel(assembly_interval, sub_interval_list, haplotagged_bam, i, hap, args):
    hap_tmp_dir = os.path.join(args.tmp_dir, "hap" + hap)
    interval_dir = os.path.join(hap_tmp_dir, "interval_" + str(i))
    threads = 1
    if not os.path.isdir(interval_dir):
        os.mkdir(interval_dir)

    for j in range(len(sub_interval_list)):
        sub_interval_dir = os.path.join(interval_dir, "sub_interval_" + str(j+1))
        if not os.path.isdir(sub_interval_dir):
            os.mkdir(sub_interval_dir)

        #######Extract interval sequence
        sub_interval = sub_interval_list[j]
        sub_interval_chr = sub_interval.split(":")[0]
        sub_interval_start = int(sub_interval.split(":")[1].split("-")[0])
        sub_interval_end = int(sub_interval.split(":")[1].split("-")[1])
        sub_interval_size = sub_interval_end - sub_interval_start + 1
        if sub_interval_size < 15000:
            sub_interval_size = 15000

        sub_interval_extend_start = 1 if sub_interval_start - 30000 < 1 else sub_interval_start - 30000
        sub_interval_extend_end = sub_interval_end + 30000
        sub_interval_extend = sub_interval_chr + ':' + str(sub_interval_extend_start) + "-" + str(sub_interval_extend_end)

        sub_interval_ref_fa = os.path.join(sub_interval_dir, "assembly_sub_interval.ref.fasta")
        sub_interval_reads_fa = os.path.join(sub_interval_dir, "assembly_sub_interval.reads.fasta")
        sub_interval_reads_bam = os.path.join(sub_interval_dir, "assembly_sub_interval.reads.bam")
        run_command("samtools faidx {} {} > {}".format(args.reference, sub_interval_extend, sub_interval_ref_fa))
        run_command("samtools view -@ {} -q {} -b {} {} | samtools fasta -F 0x100 | seqkit rmdup -j {} - -o {}".format(threads, args.min_assembly_mq, haplotagged_bam, sub_interval, threads, sub_interval_reads_fa))
        run_command("samtools view -@ {} -q {} -F 0x900 -b {} {} -o {}".format(threads, args.min_assembly_mq, haplotagged_bam, sub_interval, sub_interval_reads_bam))
        if os.path.getsize(sub_interval_reads_fa) == 0:
            #####unfinished no_trimed_reads
            continue

        #######Trim reads and Extract median reads
        sub_interval_reads_median_fa = os.path.join(sub_interval_dir, "assembly_sub_interval.reads.trim.median.fasta")
        trim_flag = 1
        try:
            run_command("mecat2map -num_threads {} -skip_overhang -db_dir {} -keep_db -task pm -outfmt m4x -out {} {} {}".format(threads, os.path.join(sub_interval_dir, "trim_pm_dir"), os.path.join(sub_interval_dir, "trim_pm_dir", "trim_pm.m4x"), sub_interval_reads_fa, sub_interval_reads_fa))
            run_command("mecat2pm4 -t {} -p 100000 -k 100 {} {} {}".format(threads, os.path.join(sub_interval_dir, "trim_pm_dir"), os.path.join(sub_interval_dir, "trim_pm4_dir"), os.path.join(sub_interval_dir, "trim_pm_dir", "trim_pm.m4x")))
            run_command("mecat2lcr -num_threads {} -out {} {} {}".format(threads, os.path.join(sub_interval_dir, "lcr.txt"), os.path.join(sub_interval_dir, "trim_pm_dir"), os.path.join(sub_interval_dir, "trim_pm4_dir")))
            run_command("mecat2splitreads -num_threads {} -out {} {} {} {}".format(threads, os.path.join(sub_interval_dir, "sr.txt"), os.path.join(sub_interval_dir, "trim_pm_dir"), os.path.join(sub_interval_dir, "trim_pm4_dir"), os.path.join(sub_interval_dir, "lcr.txt")))
            run_command("mecat2trimbases {} {} 1 {} ".format(os.path.join(sub_interval_dir, "trim_pm_dir"), os.path.join(sub_interval_dir, "sr.txt"), os.path.join(sub_interval_dir, "assembly_sub_interval.reads.trim.origin.fasta")))
            run_command("awk \'{{if(NR%2!=0) print\">\"$2; else print$0}}\' {} > {}".format(os.path.join(sub_interval_dir, "assembly_sub_interval.reads.trim.origin.fasta"), os.path.join(sub_interval_dir, "assembly_sub_interval.reads.trim.fasta")))
        except:
            trim_flag = 0

        if trim_flag:
            run_command("seqkit sort -j {} -w 0 {} | {}/pacbio_zmw_median_extract.py - {}".format(threads, os.path.join(sub_interval_dir, "assembly_sub_interval.reads.trim.fasta"), sys.path[0], sub_interval_reads_median_fa))
        else:
            run_command("seqkit sort -j {} -w 0 {} | {}/pacbio_zmw_median_extract.py - {}".format(threads, sub_interval_reads_fa, sys.path[0], sub_interval_reads_median_fa))
        if os.path.getsize(sub_interval_reads_median_fa) == 0:
            #####unfinished no_trimed_reads
            continue
        
        #######Assembly
        sub_interval_raw_contig_fa = os.path.join(sub_interval_dir, "assembly_sub_interval.raw.contig.fasta")
        if sub_interval_size <= 50000:
            run_command("minimap2 -t {} -ax map-pb --eqx --secondary=no {} {} | {}/alignment_shiftPos_toBed.py > {}".format(threads, sub_interval_ref_fa, sub_interval_reads_median_fa, sys.path[0], os.path.join(sub_interval_dir, "assembly_sub_interval.median_read_align.bed")))
            best_read_info = max_cover_read(os.path.join(sub_interval_dir, "assembly_sub_interval.median_read_align.bed"), sub_interval)
            if best_read_info[1] >= 0.95:
                run_command("seqkit grep -j {} -w 0 -r -p {} {} | awk \'{{if(NR==1)print \">ctg1\";if(NR==2)print$0}}\' > {}".format(threads, best_read_info[0], sub_interval_reads_median_fa, sub_interval_raw_contig_fa))
                #####finished representative_0.95_read
            else:
                if trim_flag:
                    run_command("wtdbg2 -t {} -x sq -g {} -S 1 --edge-min 1 --rescue-low-cov-edges --ctg-min-length {} --aln-dovetail -1 --no-read-clip -i {} -o {}".format(threads, sub_interval_size, args.min_assembly_size, sub_interval_reads_median_fa, os.path.join(sub_interval_dir, "assembly_sub_interval.contig")))
                else:
                    run_command("wtdbg2 -t {} -x sq -g {} -S 1 --edge-min 1 --rescue-low-cov-edges --ctg-min-length {} --aln-dovetail -1 -i {} -o {}".format(threads, sub_interval_size, args.min_assembly_size, sub_interval_reads_median_fa, os.path.join(sub_interval_dir, "assembly_sub_interval.contig")))
                run_command("wtpoa-cns -t {} -i {} -fo {}".format(threads, os.path.join(sub_interval_dir, "assembly_sub_interval.contig.ctg.lay.gz"), sub_interval_raw_contig_fa))
                if os.path.getsize(sub_interval_raw_contig_fa) != 0:
                    #####finished wtdbg2_assembly
                    pass
                else:
                    if best_read_info[1] >= 0.7:
                        run_command("seqkit grep -j {} -w 0 -r -p {} {} | awk \'{{if(NR==1)print \">ctg1\";if(NR==2)print$0}}\' > {}".format(threads, best_read_info[0], sub_interval_reads_median_fa, sub_interval_raw_contig_fa))
                        #####finished representative_0.7_read
                    else:
                        continue
        else:
                if trim_flag:
                    run_command("wtdbg2 -t {} -x sq -g {} -S 1 --edge-min 1 --rescue-low-cov-edges --ctg-min-length {} --aln-dovetail -1 --no-read-clip -i {} -o {}".format(threads, sub_interval_size, args.min_assembly_size, sub_interval_reads_median_fa, os.path.join(sub_interval_dir, "assembly_sub_interval.contig")))
                else:
                    run_command("wtdbg2 -t {} -x sq -g {} -S 1 --edge-min 1 --rescue-low-cov-edges --ctg-min-length {} --aln-dovetail -1 -i {} -o {}".format(threads, sub_interval_size, args.min_assembly_size, sub_interval_reads_median_fa, os.path.join(sub_interval_dir, "assembly_sub_interval.contig")))
                run_command("wtpoa-cns -t {} -i {} -fo {}".format(threads, os.path.join(sub_interval_dir, "assembly_sub_interval.contig.ctg.lay.gz"), sub_interval_raw_contig_fa))
                if os.path.getsize(sub_interval_raw_contig_fa) != 0:
                    #####finished wtdbg2_assembly
                    pass
                else:
                    #####unfinished no_raw_assembly
                    continue
        
        #######Polish
        sub_interval_racon1_contig_fa = os.path.join(sub_interval_dir, "assembly_sub_interval.racon1.contig.fasta")
        sub_interval_racon2_contig_fa = os.path.join(sub_interval_dir, "assembly_sub_interval.racon2.contig.fasta")     
        sub_interval_arrow_contig_fq = os.path.join(sub_interval_dir, "assembly_sub_interval.arrow.contig.fastq")
        sub_interval_filter_contig_fa = os.path.join(sub_interval_dir, "assembly_sub_interval.filter.contig.fasta")
        run_command("minimap2 -t {} -ax map-pb --eqx --secondary=no {} {} | samtools view -@ {} -F 0x900 -o {}".format(threads, sub_interval_raw_contig_fa, sub_interval_reads_fa, threads, os.path.join(sub_interval_dir, "assembly_sub_interval.raw.contig.sam")))
        run_command("racon -t {} -u {} {} {} > {}".format(threads, sub_interval_reads_fa, os.path.join(sub_interval_dir, "assembly_sub_interval.raw.contig.sam"), sub_interval_raw_contig_fa, sub_interval_racon1_contig_fa))
        run_command("minimap2 -t {} -ax map-pb --eqx --secondary=no {} {} | samtools view -@ {} -F 0x900 -o {}".format(threads, sub_interval_racon1_contig_fa, sub_interval_reads_fa, threads, os.path.join(sub_interval_dir, "assembly_sub_interval.racon1.contig.sam")))
        run_command("racon -t {} -u {} {} {} > {}".format(threads, sub_interval_reads_fa, os.path.join(sub_interval_dir, "assembly_sub_interval.racon1.contig.sam"), sub_interval_racon1_contig_fa, sub_interval_racon2_contig_fa))
        run_command("pbmm2 align -j {} {} {} | samtools sort -@ {} -u -O bam -o {}".format(threads, sub_interval_racon2_contig_fa, sub_interval_reads_bam, threads, os.path.join(sub_interval_dir, "assembly_sub_interval.racon2.contig.bam")))
        run_command("samtools index -@ {} {}".format(threads, os.path.join(sub_interval_dir, "assembly_sub_interval.racon2.contig.bam")))
        try: 
            run_command("gcpp -j {} -x 3 -r {} -o {} {}".format(threads, sub_interval_racon2_contig_fa, sub_interval_arrow_contig_fq, os.path.join(sub_interval_dir, "assembly_sub_interval.racon2.contig.bam")), 7200)
        except:
            run_command("gcpp -j {} -x 3 -r {} -o {} {}".format(threads*2, sub_interval_racon2_contig_fa, sub_interval_arrow_contig_fq, os.path.join(sub_interval_dir, "assembly_sub_interval.racon2.contig.bam")), 14400)
        run_command("cutadapt -j {} -q 1,1 -m 10000 {} | tr \"!\" \"+\" | seqkit seq -j {} -Q 15 | seqkit fq2fa -j {} -w 0 | awk -v j={} \'{{if(NR%2==1)print \">subinterval\"j\"_ctg\"(NR+1)/2;else print$0}}\' > {}".format(threads, sub_interval_arrow_contig_fq, threads, threads, j+1, sub_interval_filter_contig_fa))
        run_command("minimap2 -t {} -ax asm20 --eqx --secondary=no {} {} | samtools view -@ {} -h -F 0x900 | {}/alignment_shiftPos_toBed.py > {}".format(threads, sub_interval_ref_fa, sub_interval_filter_contig_fa, threads, sys.path[0], os.path.join(sub_interval_dir, "assembly_sub_interval.contig.bed")))

        if os.path.getsize(sub_interval_filter_contig_fa) == 0:
            #####unfinished no_polished_assembly
            continue

    #######interval_merge
    interval_chr = assembly_interval.split(":")[0]
    interval_start = int(assembly_interval.split(":")[1].split("-")[0])
    interval_end = int(assembly_interval.split(":")[1].split("-")[1])

    interval_raw_contig_fa = os.path.join(interval_dir, "assembly_interval.raw.contig.fasta")
    interval_merge_contig_fa = os.path.join(interval_dir, "assembly_interval.merge.contig.fasta")

    if subprocess.run("cat {}/sub_interval_*/assembly_sub_interval.filter.contig.fasta".format(interval_dir), stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell=True).returncode == 0:
        run_command("cat {}/sub_interval_*/assembly_sub_interval.filter.contig.fasta > {}".format(interval_dir, interval_raw_contig_fa))
        run_command("cat {}/sub_interval_*/assembly_sub_interval.contig.bed > {}".format(interval_dir, os.path.join(interval_dir, "assembly_interval.origin.contig.bed")))
    else:
        run_command("rm -rf {}".format(interval_dir))
        return
        
    if os.path.getsize(interval_raw_contig_fa) == 0:
        run_command("rm -rf {}".format(interval_dir))
        return


    if len(sub_interval_list) == 1:
        run_command("cat {} | seqkit seq -j {} -m 10000 -w 0 | awk -v i={} -v chr={} -v start={} -v end={} \'{{if(NR%2==1)print \">interval\"i\"_ctg\"(NR+1)/2\"_\"chr\"-\"start\"-\"end;else print$0}}\' > {}".format(interval_raw_contig_fa, threads, i, interval_chr, interval_start, interval_end, interval_merge_contig_fa))
    elif len(sub_interval_list) > 1:
        run_command("{}/phasebook_merge.py {} {} {} {}".format(os.path.join(sys.path[0], "phasebook", "scripts"), interval_raw_contig_fa, os.path.join(interval_dir, "assembly_interval.origin.contig.bed"), os.path.join(interval_dir, "contig_merge"), threads))
        run_command("cat {} | seqkit grep -j {} -r -p utg | seqkit seq -j {} -m 10000 -w 0 | awk -v i={} -v chr={} -v start={} -v end={} \'{{if(NR%2==1)print \">interval\"i\"_ctg\"(NR+1)/2\"_\"chr\"-\"start\"-\"end;else print$0}}\' > {}".format(os.path.join(interval_dir, "contig_merge", "final_contigs.rm_isolated.fa"), threads, threads, i, interval_chr, interval_start, interval_end, interval_merge_contig_fa))
    if os.path.getsize(interval_merge_contig_fa) == 0:
        run_command("rm -rf {}".format(interval_dir))
        return

    #######SV polish
    interval_reads_fa = os.path.join(interval_dir, "assembly_interval.reads.fasta")
    interval_iris_contig_fa = os.path.join(interval_dir, "assembly_interval.iris.contig.fasta")
    run_command("samtools view -@ {} -q {} -b {} {} | samtools fasta -F 0x100 | seqkit rmdup -j {} - -o {}".format(threads, args.min_assembly_mq, haplotagged_bam, assembly_interval, threads, interval_reads_fa))
    run_command("minimap2 -t {} -ax map-pb --eqx --secondary=no {} {} | samtools view -uSh -@ {} - | samtools sort -O bam -@ {} -o {}".format(threads, interval_merge_contig_fa, interval_reads_fa, threads, threads, os.path.join(interval_dir, "assembly_interval.merge.contig.bam")))
    run_command("samtools index -@ {} {}".format(threads, os.path.join(interval_dir, "assembly_interval.merge.contig.bam")))
    run_command("sniffles -t {} --input {} --reference {} --vcf {} --output-rnames --allow-overwrite".format(threads, os.path.join(interval_dir, "assembly_interval.merge.contig.bam"), interval_merge_contig_fa, os.path.join(interval_dir, "assembly_interval.sniffles.contig.vcf")))
    run_command("bcftools view --threads {} -i \"GT=\'1/1\' & SUPPORT>={} & (SVTYPE=\'INS\' | SVTYPE=\'DEL\')\" {} | bcftools view --threads {} -e \"ALT[0]==\'<INS>\'\" -o {}".format(threads, args.min_assembly_cov, os.path.join(interval_dir, "assembly_interval.sniffles.contig.vcf"), threads, os.path.join(interval_dir, "assembly_interval.sniffles.filter.contig.vcf")))
    run_command("iris genome_in={} vcf_in={} reads_in={} vcf_out={} out_dir={} threads={} --pacbio --also_deletions".format(
        interval_merge_contig_fa, os.path.join(interval_dir, "assembly_interval.sniffles.filter.contig.vcf"), os.path.join(interval_dir, "assembly_interval.merge.contig.bam"), os.path.join(interval_dir, "assembly_interval.iris.contig.vcf"), interval_dir, threads))
    run_command("bgzip -@ {} -f {} ".format(threads, os.path.join(interval_dir, "assembly_interval.iris.contig.vcf")))
    run_command("tabix -f {}".format(os.path.join(interval_dir, "assembly_interval.iris.contig.vcf.gz")))
    run_command("bcftools consensus -f {} -H 1 {} > {}".format(interval_merge_contig_fa, os.path.join(interval_dir, "assembly_interval.iris.contig.vcf.gz"), interval_iris_contig_fa))
    
    #######Contig bed
    interval_extend_start = 1 if interval_start - 30000 < 1 else interval_start - 30000
    interval_extend_end = interval_end + 30000
    interval_extend = interval_chr + ':' + str(interval_extend_start) + "-" + str(interval_extend_end)
    interval_ref_fa = os.path.join(interval_dir, "assembly_interval.ref.fasta")
    run_command("samtools faidx {} {} > {}".format(args.reference, interval_extend, interval_ref_fa))
    run_command("minimap2 -t {} -ax asm20 --eqx --secondary=no {} {} | {}/alignment_shiftPos_toBed.py > {}".format(threads, interval_ref_fa, interval_iris_contig_fa, sys.path[0], os.path.join(sub_interval_dir, "assembly_interval.contig.bed")))

    run_command("cp {} {}".format(interval_iris_contig_fa, os.path.join(hap_tmp_dir, "raw_asssembly", "assembly_interval.interval" + str(i) + ".iris.contig.fasta")))
    run_command("cp {} {}".format(os.path.join(sub_interval_dir, "assembly_interval.contig.bed"), os.path.join(hap_tmp_dir, "raw_asssembly", "assembly_interval.interval" + str(i) + ".contig.bed")))
    run_command("rm -rf {}".format(interval_dir))


def assembly(args):
    #######Assembly interval detection
    consensus_haplotag_list = os.path.join(args.out_dir, args.prefix + ".consensus.haplotag.list")
    for hap in ["1"]:
        hap_tmp_dir = os.path.join(args.tmp_dir, "hap" + hap)
        ngs_polish_dir = os.path.join(hap_tmp_dir, "ngs_polish")
        if not os.path.isdir(ngs_polish_dir):
            os.mkdir(ngs_polish_dir)

        raw_contig_fa = os.path.join(ngs_polish_dir, "assembly.raw.contig.fasta")
        raw_contig_bed = os.path.join(ngs_polish_dir, "assembly.raw.contig.bed")

        freebayes_contig_fa = os.path.join(ngs_polish_dir, "assembly.freebayes.contig.fasta")
        run_command("samtools view -@ {} -b {} -P -L {} -o {}".format(args.threads, args.ngs_bam, raw_contig_bed, os.path.join(ngs_polish_dir, "ngs.contigs_region.bam")))
        run_command("samtools view -@ {} -b -f 4 {} -o {}".format(args.threads, args.ngs_bam, os.path.join(ngs_polish_dir, "ngs.unmap.bam")))
        run_command("samtools view -@ {} -b -f 8 {} -o {}".format(args.threads, args.ngs_bam, os.path.join(ngs_polish_dir, "ngs.unmap_pair.bam")))
        run_command("samtools merge -@ {} -f -c -p -o {} {} {} {}".format(args.threads, os.path.join(ngs_polish_dir, "ngs.polish.bam"), os.path.join(ngs_polish_dir, "ngs.contigs_region.bam"), os.path.join(ngs_polish_dir, "ngs.unmap.bam"), os.path.join(ngs_polish_dir, "ngs.unmap_pair.bam")))
        run_command("samtools sort -n -@ {} {} | samtools view -@ {} -h | uniq | samtools view -@ {} -q {} -b | samtools fastq -@ {} -1 {} -2 {} -s {} > /dev/null".format(
            args.threads, os.path.join(ngs_polish_dir, "ngs.polish.bam"), args.threads, args.threads, args.min_assembly_mq, args.threads, os.path.join(ngs_polish_dir, "ngs.polish.R1.fastq"), os.path.join(ngs_polish_dir, "ngs.polish.R2.fastq"), os.path.join(ngs_polish_dir, "ngs.polish.S.fastq")))
        run_command("bwa index {}".format(raw_contig_fa))
        run_command("bwa mem -t {} -Y {} {} {} | samtools view -uSh -@ {} | samtools sort -O bam -@ {} -o {}".format(args.threads, raw_contig_fa, os.path.join(ngs_polish_dir, "ngs.polish.R1.fastq"), os.path.join(ngs_polish_dir, "ngs.polish.R2.fastq"), args.threads, args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.paired.bam")))
        run_command("bwa mem -t {} -Y {} {} | samtools view -uSh -@ {} | samtools sort -O bam -@ {} -o {}".format(args.threads, raw_contig_fa, os.path.join(ngs_polish_dir, "ngs.polish.S.fastq"), args.threads, args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.unpaired.bam")))
        run_command("samtools merge -@ {} -f -o {} {} {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.bam"), os.path.join(ngs_polish_dir, "assembly.raw.contig.unpaired.bam"), os.path.join(ngs_polish_dir, "assembly.raw.contig.paired.bam")))
        run_command("gatk MarkDuplicates -R {} -I {} -O {} -M /dev/null".format(raw_contig_fa, os.path.join(ngs_polish_dir, "assembly.raw.contig.bam"), os.path.join(ngs_polish_dir, "assembly.raw.contig.dedup.bam")))
        run_command("samtools index -@ {} {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.dedup.bam")))
        
        chroms = pysam.FastaFile(raw_contig_fa).references
        freebayes_pool = Pool(args.threads)
        for chrom in chroms:
            chrom_vcf = os.path.join(ngs_polish_dir, "assembly.freebayes." + chrom + ".vcf.gz")
            freebayes_pool.apply_async(freebayes_parallel, args = (os.path.join(ngs_polish_dir, "assembly.raw.contig.dedup.bam"), raw_contig_fa, chrom, chrom_vcf))
        freebayes_pool.close()
        freebayes_pool.join()

        run_command("find {} -name \"assembly.freebayes.*.vcf.gz\" -type f -exec ls {{}} \\; > {}".format(ngs_polish_dir, os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.list")))
        run_command("bcftools concat --threads {} --naive-force -f {} -O z -o {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.list"), os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.gz")))
        run_command("bcftools sort -O z -o {} {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.sorted.vcf.gz"), os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.gz")))
        run_command("bcftools view --threads {} -i \"GT==\'1/1\' & QUAL>1\" {} | bcftools view --threads {} -i \"FORMAT/DP>=6 & FORMAT/DP<12 & FORMAT/AD[0:0]<=0\" | {}/homopolymer_extract.py | bgzip -@ {} -c > {}".format(
            args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.sorted.vcf.gz"), args.threads, sys.path[0], args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.homopolymer.vcf.gz")))
        run_command("tabix -f {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.filter.homopolymer.vcf.gz")))
        run_command("bcftools view --threads {} -i \"GT==\'1/1\' & QUAL>1\" {} | bcftools view --threads {} -i \"FORMAT/DP>=12 & FORMAT/AD[0:0]<=0\" -o {}".format(
            args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.sorted.vcf.gz"), args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.strict.vcf.gz")))
        run_command("tabix -f {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.filter.strict.vcf.gz")))
        run_command("bcftools concat --threads {} -a {} {} -O z -o {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.strict.vcf.gz"), os.path.join(ngs_polish_dir, "assembly.freebayes.filter.homopolymer.vcf.gz"), os.path.join(ngs_polish_dir, "assembly.freebayes.filter.vcf.gz")))
        run_command("tabix -f {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.filter.vcf.gz")))
        run_command("bcftools consensus -f {} -H 1 {} > {}".format(raw_contig_fa, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.vcf.gz"), freebayes_contig_fa))
 

        #######Switch region clip
        switch_clip_dir = os.path.join(hap_tmp_dir, "switch_clip")
        if not os.path.isdir(switch_clip_dir):
            os.mkdir(switch_clip_dir)
        final_fa = os.path.join(args.out_dir, args.prefix + ".hap" + hap + ".fasta")
        run_command("minimap2 -t {} -xasm20 --cs -a {} {} | samtools sort -@ {} -O BAM -o {}".format(args.threads, args.reference, freebayes_contig_fa, args.threads, os.path.join(switch_clip_dir, "assembly.switch_clip.bam")) )
        run_command("samtools index -@ {} {}".format(args.threads, os.path.join(switch_clip_dir, "assembly.switch_clip.bam")))
        run_command("samtools view -@ {} -h {} | paftools.js sam2paf - | paftools.js call -L 10000 -f {} - | sed \"s/1\\/1/1|0/g\" > {}".format(args.threads, os.path.join(switch_clip_dir, "assembly.switch_clip.bam"), args.reference, os.path.join(switch_clip_dir, "assembly.switch_clip.vcf")))
        run_command("whatshap compare --ignore-sample-name --names truth,whatshap --switch-error-bed {} --only-snvs {} {}".format(os.path.join(switch_clip_dir, "assembly.switch_error.reference.bed"), args.vcf, os.path.join(switch_clip_dir, "assembly.switch_clip.vcf")))
        run_command("bedtools makewindows -g {} -w 50000 -s 25000 > {}".format(args.reference + '.fai', os.path.join(switch_clip_dir, "reference.window.bed")))
        run_command("bedtools intersect -a {} -b {} -wa | bedtools sort -i - | uniq -c | awk -v OFS=\'\\t\' \'{{if($1>=5)print$2,$3,$4}}\' | bedtools merge -i - > {}".format(os.path.join(switch_clip_dir, "reference.window.bed"), os.path.join(switch_clip_dir, "assembly.switch_error.reference.bed"), os.path.join(switch_clip_dir, "assembly.switch_clip.reference.bed")))
        run_command("{}/switch_region_transfer.py -c {} -b {} -o {}".format(sys.path[0], os.path.join(switch_clip_dir, "assembly.switch_clip.reference.bed"), os.path.join(switch_clip_dir, "assembly.switch_clip.bam"), os.path.join(switch_clip_dir, "assembly.switch_clip.assembly.bed")))
        run_command("bedtools sort -i {} -faidx {} | bedtools complement -i - -g {} | awk \'{{if($3-$2>=10000) print$1\":\"$2+1\"-\"$3}}\' > {}".format(os.path.join(switch_clip_dir, "assembly.switch_clip.assembly.bed"), freebayes_contig_fa + ".fai", freebayes_contig_fa + ".fai", os.path.join(switch_clip_dir, "assembly.reserve.region")))
        run_command("samtools faidx {} -r {} | seqkit seq -j {} -w 0 | awk -v sample={} -v hap=hap{} \'{{if(NR%2==1){{split($1,contig,\":\");split(contig[1],info,\"_\");split(info[3],interval,\"-\");print \">\"sample\"_\"hap\"_ctg\"(NR+1)/2\"_\"interval[1]\":\"interval[2]\"-\"interval[3]}} else print$0}}\' > {}".format(freebayes_contig_fa, os.path.join(switch_clip_dir, "assembly.reserve.region"), args.threads, args.prefix, hap, final_fa))
        run_command("rm -r {}".format(ngs_polish_dir))
        run_command("rm -r {}".format(switch_clip_dir))
        run_command("rm -r {}".format(hap_tmp_dir))

    for hap in ["2"]:
        haplotagged_read = os.path.join(args.out_dir, args.prefix + ".subreads.hap" + hap + ".reads.list")
        haplotagged_bam = os.path.join(args.out_dir, args.prefix + ".subreads.hap" + hap + ".bam")
        assembly_interval = os.path.join(args.out_dir, args.prefix + ".assembly_interval.hap" + hap + ".list")

        run_command("grep \"H{}\" {} | awk \'{{print$1}}\' > {}".format(hap, consensus_haplotag_list , haplotagged_read))
        run_command("samtools view -@ {} -b -N {} -o {} {}".format(args.threads, haplotagged_read, haplotagged_bam, args.subread_bam))
        run_command("samtools index -@ {} {}".format(args.threads, haplotagged_bam))
        run_command("rm {}".format(haplotagged_read))
        run_command("samtools depth -@ {} -Q {} -J {} | {}/assembly_interval_detect.py --min-assembly-size {} --min-assembly-cov {} --min-assembly-mq {} --min-connect-align-len {} --haplotagged-bam {} -o {}".format(
            args.threads, args.min_assembly_mq, haplotagged_bam, sys.path[0], args.min_assembly_size, args.min_assembly_cov, args.min_assembly_mq, args.min_connect_align_len, haplotagged_bam, assembly_interval))
        assembly_interval_dict = assembly_interval_split(assembly_interval, args.max_sub_interval_size, args.sub_interval_overlap_size)

        #######Interval assembly
        hap_tmp_dir = os.path.join(args.tmp_dir, "hap" + hap)
        if not os.path.isdir(hap_tmp_dir):
            os.mkdir(hap_tmp_dir)

        raw_assembly_dir = os.path.join(hap_tmp_dir, "raw_asssembly")
        if not os.path.isdir(raw_assembly_dir):
            os.mkdir(raw_assembly_dir)

        interval_assembly_pool = Pool(args.threads)
        i = 0
        for assembly_interval, sub_interval_list in assembly_interval_dict.items():
            i += 1
            interval_assembly_pool.apply_async(interval_assembly_parallel, args = (assembly_interval, sub_interval_list, haplotagged_bam, i, hap, args))
        interval_assembly_pool.close()
        interval_assembly_pool.join()


        #######Short-read polish
        ngs_polish_dir = os.path.join(hap_tmp_dir, "ngs_polish")
        if not os.path.isdir(ngs_polish_dir):
            os.mkdir(ngs_polish_dir)

        raw_contig_fa = os.path.join(ngs_polish_dir, "assembly.raw.contig.fasta")
        raw_contig_bed = os.path.join(ngs_polish_dir, "assembly.raw.contig.bed")
        run_command("> {}".format(raw_contig_fa))
        for i in range(1, len(assembly_interval_dict.keys())):
            if os.path.isfile(os.path.join(raw_assembly_dir, "assembly_interval.interval" + str(i) + ".iris.contig.fasta")):
                run_command("cat {} >> {}".format(os.path.join(raw_assembly_dir, "assembly_interval.interval" + str(i) + ".iris.contig.fasta"), raw_contig_fa))
            if os.path.isfile(os.path.join(raw_assembly_dir, "assembly_interval.interval" + str(i) + ".contig.bed")):
                run_command("cat {} >> {}".format(os.path.join(raw_assembly_dir, "assembly_interval.interval" + str(i) + ".contig.bed"), raw_contig_bed))
        run_command("rm -r {}".format(raw_assembly_dir))

        freebayes_contig_fa = os.path.join(ngs_polish_dir, "assembly.freebayes.contig.fasta")
        run_command("samtools view -@ {} -b {} -P -L {} -o {}".format(args.threads, args.ngs_bam, raw_contig_bed, os.path.join(ngs_polish_dir, "ngs.contigs_region.bam")))
        run_command("samtools view -@ {} -b -f 4 {} -o {}".format(args.threads, args.ngs_bam, os.path.join(ngs_polish_dir, "ngs.unmap.bam")))
        run_command("samtools view -@ {} -b -f 8 {} -o {}".format(args.threads, args.ngs_bam, os.path.join(ngs_polish_dir, "ngs.unmap_pair.bam")))
        run_command("samtools merge -@ {} -f -c -p -o {} {} {} {}".format(args.threads, os.path.join(ngs_polish_dir, "ngs.polish.bam"), os.path.join(ngs_polish_dir, "ngs.contigs_region.bam"), os.path.join(ngs_polish_dir, "ngs.unmap.bam"), os.path.join(ngs_polish_dir, "ngs.unmap_pair.bam")))
        run_command("samtools sort -n -@ {} {} | samtools view -@ {} -h | uniq | samtools view -@ {} -q {} -b | samtools fastq -@ {} -1 {} -2 {} -s {} > /dev/null".format(
            args.threads, os.path.join(ngs_polish_dir, "ngs.polish.bam"), args.threads, args.threads, args.min_assembly_mq, args.threads, os.path.join(ngs_polish_dir, "ngs.polish.R1.fastq"), os.path.join(ngs_polish_dir, "ngs.polish.R2.fastq"), os.path.join(ngs_polish_dir, "ngs.polish.S.fastq")))
        run_command("bwa index {}".format(raw_contig_fa))
        run_command("bwa mem -t {} -Y {} {} {} | samtools view -uSh -@ {} | samtools sort -O bam -@ {} -o {}".format(args.threads, raw_contig_fa, os.path.join(ngs_polish_dir, "ngs.polish.R1.fastq"), os.path.join(ngs_polish_dir, "ngs.polish.R2.fastq"), args.threads, args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.paired.bam")))
        run_command("bwa mem -t {} -Y {} {} | samtools view -uSh -@ {} | samtools sort -O bam -@ {} -o {}".format(args.threads, raw_contig_fa, os.path.join(ngs_polish_dir, "ngs.polish.S.fastq"), args.threads, args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.unpaired.bam")))
        run_command("samtools merge -@ {} -f -o {} {} {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.bam"), os.path.join(ngs_polish_dir, "assembly.raw.contig.unpaired.bam"), os.path.join(ngs_polish_dir, "assembly.raw.contig.paired.bam")))
        run_command("gatk MarkDuplicates -R {} -I {} -O {} -M /dev/null".format(raw_contig_fa, os.path.join(ngs_polish_dir, "assembly.raw.contig.bam"), os.path.join(ngs_polish_dir, "assembly.raw.contig.dedup.bam")))
        run_command("samtools index -@ {} {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.raw.contig.dedup.bam")))
        
        chroms = pysam.FastaFile(raw_contig_fa).references
        freebayes_pool = Pool(args.threads)
        for chrom in chroms:
            chrom_vcf = os.path.join(ngs_polish_dir, "assembly.freebayes." + chrom + ".vcf.gz")
            freebayes_pool.apply_async(freebayes_parallel, args = (os.path.join(ngs_polish_dir, "assembly.raw.contig.dedup.bam"), raw_contig_fa, chrom, chrom_vcf))
        freebayes_pool.close()
        freebayes_pool.join()

        run_command("find {} -name \"assembly.freebayes.*.vcf.gz\" -type f -exec ls {{}} \\; > {}".format(ngs_polish_dir, os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.list")))
        run_command("bcftools concat --threads {} --naive-force -f {} -O z -o {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.list"), os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.gz")))
        run_command("bcftools sort -O z -o {} {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.sorted.vcf.gz"), os.path.join(ngs_polish_dir, "assembly.freebayes.vcf.gz")))
        run_command("bcftools view --threads {} -i \"GT==\'1/1\' & QUAL>1\" {} | bcftools view --threads {} -i \"FORMAT/DP>=6 & FORMAT/DP<12 & FORMAT/AD[0:0]<=0\" | {}/homopolymer_extract.py | bgzip -@ {} -c > {}".format(
            args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.sorted.vcf.gz"), args.threads, sys.path[0], args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.homopolymer.vcf.gz")))
        run_command("tabix -f {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.filter.homopolymer.vcf.gz")))
        run_command("bcftools view --threads {} -i \"GT==\'1/1\' & QUAL>1\" {} | bcftools view --threads {} -i \"FORMAT/DP>=12 & FORMAT/AD[0:0]<=0\" -o {}".format(
            args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.sorted.vcf.gz"), args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.strict.vcf.gz")))
        run_command("tabix -f {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.filter.strict.vcf.gz")))
        run_command("bcftools concat --threads {} -a {} {} -O z -o {}".format(args.threads, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.strict.vcf.gz"), os.path.join(ngs_polish_dir, "assembly.freebayes.filter.homopolymer.vcf.gz"), os.path.join(ngs_polish_dir, "assembly.freebayes.filter.vcf.gz")))
        run_command("tabix -f {}".format(os.path.join(ngs_polish_dir, "assembly.freebayes.filter.vcf.gz")))
        run_command("bcftools consensus -f {} -H 1 {} > {}".format(raw_contig_fa, os.path.join(ngs_polish_dir, "assembly.freebayes.filter.vcf.gz"), freebayes_contig_fa))

        #######Switch region clip
        switch_clip_dir = os.path.join(hap_tmp_dir, "switch_clip")
        if not os.path.isdir(switch_clip_dir):
            os.mkdir(switch_clip_dir)
        final_fa = os.path.join(args.out_dir, args.prefix + ".hap" + hap + ".fasta")
        run_command("minimap2 -t {} -xasm20 --cs -a {} {} | samtools sort -@ {} -O BAM -o {}".format(args.threads, args.reference, freebayes_contig_fa, args.threads, os.path.join(switch_clip_dir, "assembly.switch_clip.bam")) )
        run_command("samtools index -@ {} {}".format(args.threads, os.path.join(switch_clip_dir, "assembly.switch_clip.bam")))
        run_command("samtools view -@ {} -h {} | paftools.js sam2paf - | paftools.js call -L 10000 -f {} - | sed \"s/1\\/1/1|0/g\" > {}".format(args.threads, os.path.join(switch_clip_dir, "assembly.switch_clip.bam"), args.reference, os.path.join(switch_clip_dir, "assembly.switch_clip.vcf")))
        run_command("whatshap compare --ignore-sample-name --names truth,whatshap --switch-error-bed {} --only-snvs {} {}".format(os.path.join(switch_clip_dir, "assembly.switch_error.reference.bed"), args.vcf, os.path.join(switch_clip_dir, "assembly.switch_clip.vcf")))
        run_command("bedtools makewindows -g {} -w 50000 -s 25000 > {}".format(args.reference + '.fai', os.path.join(switch_clip_dir, "reference.window.bed")))
        run_command("bedtools intersect -a {} -b {} -wa | bedtools sort -i - | uniq -c | awk -v OFS=\'\\t\' \'{{if($1>=5)print$2,$3,$4}}\' | bedtools merge -i - > {}".format(os.path.join(switch_clip_dir, "reference.window.bed"), os.path.join(switch_clip_dir, "assembly.switch_error.reference.bed"), os.path.join(switch_clip_dir, "assembly.switch_clip.reference.bed")))
        run_command("{}/switch_region_transfer.py -c {} -b {} -o {}".format(sys.path[0], os.path.join(switch_clip_dir, "assembly.switch_clip.reference.bed"), os.path.join(switch_clip_dir, "assembly.switch_clip.bam"), os.path.join(switch_clip_dir, "assembly.switch_clip.assembly.bed")))
        run_command("bedtools sort -i {} -faidx {} | bedtools complement -i - -g {} | awk \'{{if($3-$2>=10000) print$1\":\"$2+1\"-\"$3}}\' > {}".format(os.path.join(switch_clip_dir, "assembly.switch_clip.assembly.bed"), freebayes_contig_fa + ".fai", freebayes_contig_fa + ".fai", os.path.join(switch_clip_dir, "assembly.reserve.region")))
        run_command("samtools faidx {} -r {} | seqkit seq -j {} -w 0 | awk -v sample={} -v hap=hap{} \'{{if(NR%2==1){{split($1,contig,\":\");split(contig[1],info,\"_\");split(info[3],interval,\"-\");print \">ctg\"(NR+1)/2\"_\"interval[1]\"-\"interval[2]\"-\"interval[3]}} else print$0}}\' > {}".format(freebayes_contig_fa, os.path.join(switch_clip_dir, "assembly.reserve.region"), args.threads, args.prefix, hap, final_fa))
        run_command("rm -r {}".format(ngs_polish_dir))
        run_command("rm -r {}".format(switch_clip_dir))
        run_command("rm -r {}".format(hap_tmp_dir))




def main():
    args = parse_args()
    check_software()
    check_file(args.reference)
    check_file(args.subread_bam)
    check_file(args.hifi_bam)
    check_file(args.vcf)
    check_file(args.par_bed)
    if args.gender != "male" and args.gender != "female":
        logger.error("Invalid gender: " + args.gender)
        sys.exit(1)

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    if not os.path.isdir(args.tmp_dir):
        os.mkdir(args.tmp_dir)

    #haplotag(args)
    assembly(args)

if __name__ == "__main__":
    main()

