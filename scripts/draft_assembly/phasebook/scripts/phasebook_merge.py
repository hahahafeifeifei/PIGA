#!/usr/bin/env python3

# Phasebook merge script
#
# Utlize the assemble_supereads module of phasebook to assemble the interval contigs.
# 
# Date: 2022/12/10

import sys
import os
from assembly import assemble_supereads

asm_supereads_fa = sys.argv[1]
asm_supereads_bed = sys.argv[2]
asm_supereads_dir = sys.argv[3]
threads = int(sys.argv[4])

if not os.path.isdir(asm_supereads_dir):
    os.mkdir(asm_supereads_dir)

min_read_len=10000
sp_min_ovlplen=15000
sp_min_identity=0.97
sp_oh=5000
sp_ohratio=0.3
max_tip_len=1000
ctg_asm='naive'
max_het_snps=0
min_allele_cov=4
platform='pb'
rm_tmp=True
super_ovlp_fast=False

assemble_supereads(asm_supereads_fa, asm_supereads_dir, asm_supereads_bed, threads, min_read_len, sp_min_ovlplen, sp_min_identity, sp_oh, sp_ohratio, max_tip_len, ctg_asm, max_het_snps, min_allele_cov, platform, rm_tmp, super_ovlp_fast)
