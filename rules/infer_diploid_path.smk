
rule all_infer_diploid_path:
    input:
        expand("c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.clip.fasta", sample=config['samples']),
        expand("c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.clip.fasta", sample=config['samples'])

rule kmc_prepare_sample_kmer:
    input:
        ngs_R1_fq = lambda wildcards: get_sr_input_fastqs(wildcards)[0],
        ngs_R2_fq = lambda wildcards: get_sr_input_fastqs(wildcards)[1],
        hifi_fq = get_hifi_input_fastqs
    output:
        kmer_list = "c8_diploid_path_infer/result/{sample}/{sample}.kmer.list",
        kff = temp("c8_diploid_path_infer/result/{sample}/{sample}.kmer.kff")
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        > {output.kmer_list}
        echo {input.ngs_R1_fq} >> {output.kmer_list}
        echo {input.ngs_R2_fq} >> {output.kmer_list}
        echo {input.hifi_fq} >> {output.kmer_list}
        
        mkdir c8_diploid_path_infer/result/{wildcards.sample}/tmp
        kmc -k29 -m28 -o kff -t {threads} -hp @{output.kmer_list} {output.kff} c8_diploid_path_infer/result/{wildcards.sample}/tmp
        rm -rf c8_diploid_path_infer/result/{wildcards.sample}/tmp
        """
        
rule get_sample_haplotype_path:
    input:
        hapl = get_hapl_input,
        gbz = get_gbz_input,
        sample_kff = "c8_diploid_path_infer/result/{sample}/{sample}.kmer.kff"
    output:
        sample_haplotype_gbz = temp("c8_diploid_path_infer/result/{sample}/{sample}.haplotypes.gbz"),
        sample_haplotype_gfa = temp("c8_diploid_path_infer/result/{sample}/{sample}.haplotypes.gfa")
    threads: 4
    resources:
        #450G
        mem_mb = 450000
    shell:
        """
        vg haplotypes -t {threads} -i {input.hapl} -k {input.sample_kff} -g {output.sample_haplotype_gbz} {input.gbz} --num-haplotypes 17
        
        vg convert -f {output.sample_haplotype_gbz} > {output.sample_haplotype_gfa}
        
        """
#TODO: how to get the chr.node_range.list?
rule sample_haplotype_mod:
    input:
        range = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/chr.node_range.list",
        sample_haplotype_gfa = "c8_diploid_path_infer/result/{sample}/{sample}.haplotypes.gfa"
    output:
        sample_haplotype_mod_gfa = temp("c8_diploid_path_infer/result/{sample}/{sample}.haplotypes.mod.gfa")
    threads: 1
    resources:
        #20G
        mem_mb=20000
    shell:
        """
        python3 ~/software/script/graph-genotyping/gfa_haplotype_mod.py \
        {input.sample_haplotype_gfa} \
        {input.range} \
        {output.sample_haplotype_mod_gfa}
        """
        
rule form_sample_merge_gfa_auto_X:
    input:
        gfa = get_gfa_input,
        variant_path = get_variant_path_input,
        sample_haplotype_mod_gfa = "c8_diploid_path_infer/result/{sample}/{sample}.haplotypes.mod.gfa"
    output:
        sample_merge_gfa = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.gfa")
    wildcard_constraints:
        chr = "chr([1-9]|1[0-9]|2[0-2]|X)"
    shell:
        """
        > {output.sample_merge_gfa}
        grep "^S\\|^L\\|{wildcards.sample}" {input.gfa} >> {output.sample_merge_gfa}
        
        grep "{wildcards.sample}" {input.variant_path} >> {output.sample_merge_gfa}
        
        grep -P "{wildcards.chr}\\t" {input.sample_haplotype_mod_gfa} >> {output.sample_merge_gfa}
        """

rule form_sample_merge_gfa_Y_M:
    input:
        gfa = get_gfa_input,
        variant_path = get_variant_path_input,
        sample_haplotype_mod_gfa = "c8_diploid_path_infer/result/{sample}/{sample}.haplotypes.mod.gfa"
    output:
        sample_merge_gfa = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.gfa")
    wildcard_constraints:
        chr = "chr(Y|M)"  
    shell:
        """
        > {output.sample_merge_gfa}
        grep "^S\\|^L\\|{wildcards.sample}" {input.gfa} >> {output.sample_merge_gfa}
        
        grep -P "{wildcards.chr}\\t" {input.sample_haplotype_mod_gfa} >> {output.sample_merge_gfa}
        """

#TODO: copy the script into github.
rule sample_gfa_deconstruct:
    input:
        sample_merge_gfa = "c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.gfa"
    output:
        rmac0_gfa = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.rmac0.gfa"),
        ref_fasta = "c8_diploid_path_infer/result/{sample}/{sample}.{chr}.ref.fasta",
        unchop_vg = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.unchop.vg"),
        unchop_gfa = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.unchop.gfa"),
        unchop_gfaffix_gfa = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.unchop.gfaffix.gfa"),
        unchop_gfaffix_info = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.merge.unchop.gfaffix.info"),
        vcf = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.vcf")
    wildcard_constraints:
        chr = "chr([1-9]|1[0-9]|2[0-2]|X|Y|M)"
    threads: 4
    resources:
        #85G
        mem_mb = 85000
    shell:
        """
        python3 scripts/graph-simplification/gfa_remove_ac0_reverse.py {input.sample_merge_gfa} {output.rmac0_gfa}

        python3 scripts/graph-genotyping/gfa_fa.py {output.rmac0_gfa} recombination_ref 0 {output.ref_fasta}
        
        vg mod -u {output.rmac0_gfa} > {output.unchop_vg}

        vg view {output.unchop_vg} > {output.unchop_gfa}
        
        gfaffix {output.unchop_gfa} -o {output.unchop_gfaffix_gfa} > {output.unchop_gfaffix_info}

        vg deconstruct -t {threads} -P recombination_ref -a -e -C {output.unchop_gfaffix_gfa} > {output.vcf}

        """

rule ref_fasta_concat:
    input:
        ref_fastas = expand("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.ref.fasta", chr=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"], allow_missing=True)
    output:
        concat_ref_fasta = "c8_diploid_path_infer/result/{sample}/{sample}.ref.fasta"
    params:
        sex = get_sex
    shell:
        """
        if [ {params.sex} -eq "female" ];
        then \
                cat c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.chr{{{{1..22}},X,M}}.ref.fasta > {output.concat_ref_fasta}
        else \
                cat c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.chr{{{{1..22}},X,Y,M}}.ref.fasta > {output.concat_ref_fasta}
        fi        
        rm c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.chr{{{{1..22}},X,Y,M}}.ref.fasta
        """
        
rule chr_vcf_norm:
    input:
        vcf = "c8_diploid_path_infer/result/{sample}/{sample}.{chr}.vcf"
    output:
        norm_vcf = temp("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.norm.vcf.gz")
    params:
        sex = get_sex
    threads: 4
    shell:
        """
        python3 ~/software/script/graph-genotyping/pangenie_norm.py {input.vcf} {wildcards.sample},{wildcards.sample}.snv,recombination {wildcards.sample} {params.sex} > c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.{wildcards.chr}.norm.vcf
        bgzip -f -@ {threads} c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.{wildcards.chr}.norm.vcf
        tabix -f {output.norm_vcf}
        """
        
    
rule chr_vcf_concat:
    input:
        norm_chr_vcfs = expand("c8_diploid_path_infer/result/{sample}/{sample}.{chr}.norm.vcf.gz", chr=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"], allow_missing=True)
    output:
        pan_vcf = temp("c8_diploid_path_infer/result/{sample}/{sample}.pan.vcf.gz")
    threads: 4
    shell:
        """
        bcftools concat --threads {threads} {input.norm_chr_vcfs} -o {output.pan_vcf}
        tabix -f {output.pan_vcf}
        """

rule pan_vcf_vcfbub:
    input:
        pan_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan.vcf.gz"
    output:
        vcfbub_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan.vcfbub.vcf"
    shell:
        """
        vcfbub --input {input.pan_vcf} -l 0 -a 50000 | bcftools view -a - | bcftools view -e "ALT=='.'" -o {output.vcfbub_vcf}

        """
        
rule kmer_fq2fa:
    input:
        kmer_list = "c8_diploid_path_infer/result/{sample}/{sample}.kmer.list"
    output:
        kmer_fa = "c8_diploid_path_infer/result/{sample}/{sample}.kmer.fa"
    threads: 4
    shell:
        """
        cat {input.kmer_list} | while read line;
        do \
            seqkit fq2fa -j {threads} $line;
        done > {output.kmer_fa}
        """  

### genotype

rule vcfbub_vcf_filter:
    input:
        vcfbub_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan.vcfbub.vcf"
    output:
        vcfbub_fiter_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan.vcfbub.fiter.vcf",
        vcfbub_bed = temp("c8_diploid_path_infer/result/{sample}/{sample}.pan.vcfbub.bed"),
        vcfbub_vcf_to_filter = temp("c8_diploid_path_infer/result/{sample}/{sample}.pan.vcfbub.vcf.to_filter")
    shell:
        """
        grep -v "^#" {wildcards.sample}.pan.vcfbub.vcf | awk '{{start=$2-1;len=length($4);end=start+len;print $1"\\t"start"\\t"end"\\t"$3}}' > {wildcards.sample}.pan.vcfbub.vcf.bed
        bedtools intersect -a {wildcards.sample}.pan.vcfbub.vcf.bed -b {wildcards.sample}.pan.vcfbub.vcf.bed -wa | uniq -c | awk '{{if($1>=10)print $5}}' > {wildcards.sample}.pan.vcfbub.vcf.to_filter
        grep -v -f {wildcards.sample}.pan.vcfbub.vcf.filter {wildcards.sample}.pan.vcfbub.vcf > {wildcards.sample}.pan.vcfbub.filter.vcf 
        
        
        """  
rule vcfbub_vcf_genotype:
    input:
        vcfbub_fiter_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan.vcfbub.fiter.vcf",
        concat_ref_fasta = "c8_diploid_path_infer/result/{sample}/{sample}.ref.fasta",
        kmer_fa = "c8_diploid_path_infer/result/{sample}/{sample}.kmer.fa"
    output:
        genotype_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan_genotyping.vcf.gz",
        phasing_vcf = temp("c8_diploid_path_infer/result/{sample}/{sample}.pan_phasing.vcf.gz")
    threads: 4
    shell:
        """
        ~/software/pangenie/build/src/PanGenie-index -v {input.vcfbub_fiter_vcf} -r {input.concat_ref_fasta} -t {threads} -o c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pangenie
        ~/software/pangenie/build/src/PanGenie -j {threads} -t {threads} -i {input.kmer_fa} -f c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pangenie -o c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan
        ~/software/pangenie/build/src/PanGenie -j {threads} -t {threads} -i {input.kmer_fa} -f c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pangenie -p -o c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan
        
        bgzip -f -@ {threads} c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_genotyping.vcf
        tabix -f {output.genotype_vcf}
        
        bgzip -f -@ {threads} c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_phasing.vcf
        tabix -f {output.phasing_vcf}
        
        rm c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pangenie* c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_histogram.histo
        """

rule pangenie_phase_fix:
    input:
        phasing_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan_phasing.vcf.gz"
    output:
        phasing_fix_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan_phasing.fix.vcf.gz"
    threads: 10
    resources:
        #80G
        mem_mb = 80000
    shell:
        """
        python3 ~/software/script/graph-simplification/pangenie_phase_fix.py {input.phasing_vcf} c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_phasing.fix.vcf
        
        bgzip -f -@ {threads} c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_phasing.fix.vcf
        tabix -f {output.phasing_fix_vcf}
        """  
        
rule genotype_refine:
    input:
        vcfbub_fiter_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan.vcfbub.fiter.vcf",
        genotype_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan_genotyping.vcf.gz"
    output:
        genotype_final_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.genotype.vcf.gz"
    shell:
        """
        bcftools merge -m all --force-samples c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_genotyping.vcf.gz c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan.vcfbub.vcf | python3 ~/software/script/graph-genotyping/graph_genotype_refine.py - c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_genotyping.refine.vcf sample {wildcards.sample}.snv {wildcards.sample}
        
        bcftools view -a c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan_genotyping.refine.vcf | bcftools view -i "ALT!='.'" -o {output.genotype_final_vcf}
        tabix -f {output.genotype_final_vcf}
        """  

#TODO: script
rule genotype_phase:
    input:
        genotype_final_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.genotype.vcf.gz",
        phasing_fix_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pan_phasing.fix.vcf.gz",
    output:
        pangenie_phase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pangenie.phase.vcf.gz",
        assembly_phase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.assembly.phase.vcf.gz",
        snv_phase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.snv.phase.vcf.gz"
    threads: 10
    resources:
        #80G
        mem_mb = 80000
    shell:
        """
        bcftools merge -m all --force-samples \
            {input.genotype_final_vcf} \
            {input.phasing_fix_vcf} | \
            python3 scripts/graph-genotyping/genotype_phase.py - c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pangenie.phase.vcf
        
        bgzip -@ {threads} -f {wildcards.sample}.pangenie.phase.vcf
        tabix -f {output.pangenie_phase_vcf}
        
        bcftools merge -m all --force-samples \
            {input.genotype_final_vcf} \
            c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan.vcfbub.vcf.gz | \
            bcftools view -s sample,{wildcards.sample} | \
            python3 scripts/graph-genotyping/genotype_phase.py - c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.assembly.phase.vcf
        
        bgzip -@ {threads} -f c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.assembly.phase.vcf
        tabix -f {output.assembly_phase_vcf}
        
        bcftools merge -m all --force-samples \
            {input.genotype_final_vcf} \
            c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.pan.vcfbub.vcf.gz | \
            bcftools view -s sample,{wildcards.sample}.snv | \
            python3 scripts/graph-genotyping/genotype_phase.py - c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.snv.phase.vcf
        
        bgzip -@ {threads} -f c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.snv.phase.vcf
        tabix -f {output.snv_phase_vcf}
        """  

rule sample_ref_realign:
    input:
        concat_ref_fasta = "c8_diploid_path_infer/result/{sample}/{sample}.ref.fasta",
        kmer_list = "c8_diploid_path_infer/result/{sample}/{sample}.kmer.list"
    output:
        zmw_bam = temp("c8_diploid_path_infer/result/{sample}/{sample}.ref.zmw.bam")
    threads: 10
    resources:
        #80G
        mem_mb = 80000
    shell:
        """
        minimap2 -t {threads} \
            -ax map-pb \
            -Y -L --eqx --cs \
            {input.concat_ref_fasta} \
            $(grep CLR {input.kmer_list} | sed 's/Q20.success/merge/g')  | samtools view -@ {threads} -hb - | samtools sort -@ {threads} -o {output.zmw_bam}
        
        samtools index -@ {threads} {output.zmw_bam}
        """  
rule hiphase_phase:
    input:
        concat_ref_fasta = "c8_diploid_path_infer/result/{sample}/{sample}.ref.fasta",
        zmw_bam = "c8_diploid_path_infer/result/{sample}/{sample}.ref.zmw.bam",
        genotype_final_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.genotype.vcf.gz"
    output:
        hiphase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hiphase.phase.vcf.gz"
    threads: 10
    resources:
        #80G
        mem_mb = 80000
    shell:
        """
        hiphase -t {threads} \
            --ignore-read-groups \
            --reference {input.concat_ref_fasta} \
            --vcf {input.genotype_final_vcf} \
            --bam {input.zmw_bam} \
            --output-vcf {output.hiphase_vcf} \
            --global-realignment-cputime 600
        
        tabix -f {output.hiphase_vcf}
        """
#TODO:why multiple attempts?        
rule margin_phase:
    input:
        concat_ref_fasta = "c8_diploid_path_infer/result/{sample}/{sample}.ref.fasta",
        zmw_bam = "c8_diploid_path_infer/result/{sample}/{sample}.ref.zmw.bam",
        genotype_final_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.genotype.vcf.gz"
    output:
        chr_margin_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.margin.{chr}.phased.vcf.gz"
    shell:
        """
        margin phase \
            {input.zmw_bam} \
            {input.concat_ref_fasta} \
            {input.genotype_final_vcf} \
            ~/software/margin/params/phase/allParams.phase_vcf.ont.sv.json \
            -M -t 1 \
            -o c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.margin.{wildcards.chr} \
            -r {wildcards.chr}
        
        bgzip -f c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.margin.{wildcards.chr}.phased.vcf
        tabix -f {output.chr_margin_vcf}
        
        """  
rule concat_margin_phased_vcf:
    input:
        chr_margin_vcfs = expand(
            "c8_diploid_path_infer/result/{sample}/{sample}.margin.{chr}.phased.vcf.gz",
            chr=get_sex_specific_chr_list,
            allow_missing=True
        )
    output:
        concat_margin_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.margin.phase.vcf.gz"
    threads: 10
    shell:
        """
        bcftools concat --threads {threads} {input.chr_margin_vcfs} -o {output.concat_margin_vcf}
        tabix {output.concat_margin_vcf}
        
        rm c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.margin.chr*.chunks.csv c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.margin.chr*.phaseset.bed
        rm c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.margin.chr*.phased.vcf.gz c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.margin.chr*.phased.vcf.gz.tbi

        """  
        
rule integrate_phase:
    input:
        genotype_final_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.genotype.vcf.gz",
        pangenie_phase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.pangenie.phase.vcf.gz",
        assembly_phase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.assembly.phase.vcf.gz",
        snv_phase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.snv.phase.vcf.gz",
        hiphase_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hiphase.phase.vcf.gz",
        concat_margin_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.margin.phase.vcf.gz"
    output:
        integrate_phase_vcf = temp("c8_diploid_path_infer/result/{sample}/{sample}.integrate.phase.vcf.gz"),
        integrate_phase_filter_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.integrate.phase.filter.vcf.gz"
    shell:
        """
        bcftools merge -m all --force-samples \
            {input.genotype_final_vcf} \
            {input.snv_phase_vcf} \
            {input.assembly_phase_vcf} \
            {input.pangenie_phase_vcf} \
            {input.hiphase_vcf} \
            {input.concat_margin_vcf} | \
            python3 ~/software/script/graph-genotyping/graph_phase_intergrate.py - {output.integrate_phase_vcf} sample 2:sample,3:sample,4:sample,5:sample,6:sample > c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.integrate.phase.log
        
        bcftools view -a {output.integrate_phase_vcf} | bcftools view -i "ALT!='.'" -o {output.integrate_phase_filter_vcf}
        tabix {output.integrate_phase_filter_vcf}

        """  
# consensus       
rule complete_assembly:
    input:
        concat_ref_fasta = "c8_diploid_path_infer/result/{sample}/{sample}.ref.fasta",
        integrate_phase_filter_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.integrate.phase.filter.vcf.gz"
    output:
        complete_assembly_hap1_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.fasta"
    params:
        sex = get_sex
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        if [ {params.sex} -eq "female" ]; then
            > {output.complete_assembly_hap1_fa}
            samtools faidx {input.concat_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
            samtools faidx {input.concat_ref_fasta} chrX | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
                
            > {output.complete_assembly_hap2_fa}
            samtools faidx {input.concat_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 2 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
            samtools faidx {input.concat_ref_fasta} chr{{X,M}} | \
                bcftools consensus -H 2 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
                
        else 
            > {output.complete_assembly_hap1_fa}
            samtools faidx {input.concat_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
            samtools faidx {input.concat_ref_fasta} chrY | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
                
            > {output.complete_assembly_hap2_fa}
            samtools faidx {input.concat_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 2 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
            samtools faidx {input.concat_ref_fasta} chr{{X,M}} | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
        fi
        """  
        

        
rule merge_complete_assembly:
    input:
        complete_assembly_hap1_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.fasta"
    output:
        merge_complete_assembly_fa = temp("c8_diploid_path_infer/result/{sample}/{sample}.merge.complete_assembly.fasta")
    shell:
        """
        > {output.merge_complete_assembly_fa}
        sed "s/chr/H1_chr/g" {input.complete_assembly_hap1_fa} >> {output.merge_complete_assembly_fa}                                                                       
        sed "s/chr/H2_chr/g" {input.complete_assembly_hap2_fa} >> {output.merge_complete_assembly_fa}
        """
### coorect by phase assembly
      
rule secphase_correct_bam:
    input:
        merge_complete_assembly_fa = "c8_diploid_path_infer/result/{sample}/{sample}.merge.complete_assembly.fasta",
        phase_assembly_hap1 = get_hap1_adaptor_masked_fa_input,
        phase_assembly_hap2 = get_hap2_adaptor_masked_fa_input
    output:
        complete_assembly_phase_assembly_bam = temp("c8_diploid_path_infer/result/{sample}/{sample}.merge.complete_assembly.phase_assembly.bam"),
        secphase_log = temp("c8_diploid_path_infer/result/{sample}/secphase_out_dir/secphase.out.log"),
        complete_assembly_correct_bam = temp("c8_diploid_path_infer/result/{sample}/{sample}.merge.complete_assembly.phase_assembly.correct.bam")
    threads: 4
    shell:
        """
        minimap2 -t {threads} \
            -I 8G \
            -ax asm20 \
            -Y -L --eqx --cs \
            {input.merge_complete_assembly_fa}
            {input.phase_assembly_hap1} \
            {input.phase_assembly_hap2} | samtools sort -@ {threads} -n | samtools view -@ {threads} -hb \
            > {output.complete_assembly_phase_assembly_bam}
        
        secphase -@ {threads} --hifi -i {output.complete_assembly_phase_assembly_bam} -f {input.merge_complete_assembly_fa} -o c8_diploid_path_infer/result/{wildcards.sample}/secphase_out_dir
        
        correct_bam --threads {threads} -i {output.complete_assembly_phase_assembly_bam} -P {output.secphase_log} -o {output.complete_assembly_correct_bam} --primaryOnly -m 0 -a 0
        """
        
rule get_phase_assembly_correct_fasta:
    input:
        complete_assembly_correct_bam = "c8_diploid_path_infer/result/{sample}/{sample}.merge.complete_assembly.phase_assembly.correct.bam",
    output:
        phase_assembly_correct_hap1 = "c8_diploid_path_infer/result/{sample}/{sample}.phase_assembly.H1.fasta",
        phase_assembly_correct_hap2 = "c8_diploid_path_infer/result/{sample}/{sample}.phase_assembly.H2.fasta"
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        samtools view -@ {threads} {input.complete_assembly_correct_bam} | grep H1 | awk '{{print ">"$1;print$10}}' | seqkit rmdup -o {output.phase_assembly_correct_hap1}
        
        samtools view -@ {threads} {input.complete_assembly_correct_bam} | grep H2 | awk '{{print ">"$1;print$10}}' | seqkit rmdup -o {output.phase_assembly_correct_hap2}
        """
        
rule get_phase_assembly_vcf:
    input:
        complete_assembly_hap1_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.fasta",
        phase_assembly_correct_hap1 = "c8_diploid_path_infer/result/{sample}/{sample}.phase_assembly.H1.fasta",
        phase_assembly_correct_hap2 = "c8_diploid_path_infer/result/{sample}/{sample}.phase_assembly.H2.fasta"
    output:
        complete_assembly_phase_assembly_hap1_vcf = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.vcf.gz"),
        complete_assembly_phase_assembly_hap2_vcf = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.vcf.gz")
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        minimap2 -t {threads} \
            -c --cs -x asm20 -B 2 -E 3,1 -O 6,100 \
            {input.complete_assembly_hap1_fa} \
            {input.phase_assembly_correct_hap1} | sort -k6,6 -k8,8n | paftools.js call \
            -f {input.complete_assembly_hap1_fa} -l 5000 -L 5000 - \
            > c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.hap1.complete_assembly.phase_assembly.vcf
        bgzip -@ {threads} -f c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.hap1.complete_assembly.phase_assembly.vcf
        tabix -f {output.complete_assembly_phase_assembly_hap1_vcf}
        
        minimap2 -t {threads} \
            -c --cs -x asm20 -B 2 -E 3,1 -O 6,100 \
            {input.complete_assembly_hap2_fa} \
            {input.phase_assembly_correct_hap2} | sort -k6,6 -k8,8n | paftools.js call \
            -f {input.complete_assembly_hap2_fa} -l 5000 -L 5000 - \
            > c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.hap2.complete_assembly.phase_assembly.vcf
        bgzip -@ {threads} -f c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.hap2.complete_assembly.phase_assembly.vcf
        tabix -f {output.complete_assembly_phase_assembly_hap2_vcf}
        """
        
rule get_phase_assembly_sv_vcf:
    input:
        complete_assembly_phase_assembly_hap1_vcf = "{sample}.hap1.complete_assembly.phase_assembly.vcf.gz",
        complete_assembly_phase_assembly_hap2_vcf = "{sample}.hap2.complete_assembly.phase_assembly.vcf.gz"
    output:
        complete_assembly_phase_assembly_hap1_sv = temp("{sample}.hap1.complete_assembly.phase_assembly.sv.vcf"),
        complete_assembly_phase_assembly_hap2_sv = temp("{sample}.hap2.complete_assembly.phase_assembly.sv.vcf")
    shell:
        """
        zcat {input.complete_assembly_phase_assembly_hap1_vcf} | awk '{{if(substr($0,1,1)=="#")print$0;else{{if(length($5)-length($4)>=50 || length($4)-length($5)>=50)print$0}} }}' | \
            truvari anno svinfo > {output.complete_assembly_phase_assembly_hap1_sv}
        
        zcat {input.complete_assembly_phase_assembly_hap2_vcf} | awk '{{if(substr($0,1,1)=="#")print$0;else{{if(length($5)-length($4)>=50 || length($4)-length($5)>=50)print$0}} }}' | \
            truvari anno svinfo > {output.complete_assembly_phase_assembly_hap2_sv}
        """
### correct by raw reads        
rule get_complete_assembly_zmw_reads_bam:
    input:
        complete_assembly_hap1_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.fasta",
        kmer_list = "c8_diploid_path_infer/result/{sample}/{sample}.kmer.list"
    output:
        complete_assembly_zmw_reads_hap1_bam = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.zmw.bam"),
        complete_assembly_zmw_reads_hap2_bam = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.zmw.bam")
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        minimap2 -t {threads} \
            -ax map-pb \
            -Y -L --eqx --cs \
            {input.complete_assembly_hap1_fa} \
            $(grep CLR {input.kmer_list} | sed 's/Q20.success/merge/g') | samtools view -@ {threads} -hb - | \
            samtools sort -@ {threads} -o {output.complete_assembly_zmw_reads_hap1_bam}
        
        samtools index -@ {threads} {output.complete_assembly_zmw_reads_hap1_bam}
         
        minimap2 -t {threads} \
            -ax map-pb \
            -Y -L --eqx --cs \
            {input.complete_assembly_hap2_fa} \
            $(grep CLR {input.kmer_list} | sed 's/Q20.success/merge/g') | samtools view -@ {threads} -hb - | \
            samtools sort -@ {threads} -o {output.complete_assembly_zmw_reads_hap2_bam}
        
        samtools index -@ {threads} {output.complete_assembly_zmw_reads_hap2_bam}
        """
        
        
        
rule cuteSV_call_sv:
    input:
        complete_assembly_hap1_zmw_reads_bam = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.zmw.bam",
        complete_assembly_hap2_zmw_reads_bam = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.zmw.bam",
        complete_assembly_hap1_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.fasta",
        complete_assembly_phase_assembly_hap1_sv = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.sv.vcf",
        complete_assembly_phase_assembly_hap2_sv = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.sv.vcf"
    output:
        complete_assembly_hap1_cuteSV = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.cutesv.vcf"),
        complete_assembly_hap2_cuteSV = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.cutesv.vcf"),
        complete_assembly_hap1_true_cuteSV = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.cutesv.true.vcf"),
        complete_assembly_hap2_true_cuteSV = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.cutesv.true.vcf")
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        mkdir c8_diploid_path_infer/result/{wildcards.sample}/cutesv.tmp
        
        cuteSV \
            {input.complete_assembly_hap1_zmw_reads_bam} \
            {input.complete_assembly_hap1_fa} \
            {output.complete_assembly_hap1_cuteSV} \
            c8_diploid_path_infer/result/{wildcards.sample}/cutesv.tmp \
            -Ivcf {input.complete_assembly_phase_assembly_hap1_sv} \
            -t {threads} \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 200 \
            --diff_ratio_merging_DEL 0.5 \
            -q 10 -L -1
            
        rm -r c8_diploid_path_infer/result/{wildcards.sample}/cutesv.tmp
        bcftools view -i "GT!='RR'" {output.complete_assembly_hap1_cuteSV} -o {output.complete_assembly_hap1_true_cuteSV}
        
        mkdir c8_diploid_path_infer/result/{wildcards.sample}/cutesv.tmp
        
        cuteSV \
            {input.complete_assembly_hap2_zmw_reads_bam} \
            {input.complete_assembly_hap2_fa} \
            {output.complete_assembly_hap2_cuteSV} \
            c8_diploid_path_infer/result/{wildcards.sample}/cutesv.tmp \
            -Ivcf {input.complete_assembly_phase_assembly_hap2_sv} \
            -t {threads} \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 200 \
            --diff_ratio_merging_DEL 0.5 \
            -q 10 -L -1
            
        rm -r c8_diploid_path_infer/result/{wildcards.sample}/cutesv.tmp
        bcftools view -i "GT!='RR'" {output.complete_assembly_hap2_cuteSV} -o {output.complete_assembly_hap2_true_cuteSV}

        """
rule complete_assembly_get_polish_region:
    input:
        complete_assembly_hap1_true_cuteSV = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.cutesv.true.vcf",
        complete_assembly_hap2_true_cuteSV = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.cutesv.true.vcf"
    output:
        complete_assembly_hap1_polish_region_bed = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.polish_region.bed"),
        complete_assembly_hap2_polish_region_bed = temp("c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.polish_region.bed")
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        bcftools query -f '%CHROM\\t%POS\\t%END\\t%SVLEN\\n' {input.complete_assembly_hap1_true_cuteSV} | awk '{{if($4<=0) svlen=-1*$4;else svlen=$4;if(($2-3*svlen)<=0)print $1"\\t0\\t"$2+3*svlen;else print $1"\\t"$2-3*svlen"\\t"$2+3*svlen}}' | bedtools sort -i | bedtools merge -i - > {output.complete_assembly_hap1_polish_region_bed}
        
        bcftools query -f '%CHROM\\t%POS\\t%END\\t%SVLEN\\n' {input.complete_assembly_hap2_true_cuteSV} | awk '{{if($4<=0) svlen=-1*$4;else svlen=$4;if(($2-3*svlen)<=0)print $1"\\t0\\t"$2+3*svlen;else print $1"\\t"$2-3*svlen"\\t"$2+3*svlen}}' | bedtools sort -i | bedtools merge -i - > {output.complete_assembly_hap2_polish_region_bed}
        

        """
rule complete_assembly_get_polish_region_vcf:
    input:
        complete_assembly_hap1_polish_region_bed = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.polish_region.bed",
        complete_assembly_hap2_polish_region_bed = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.polish_region.bed",
        complete_assembly_phase_assembly_hap1_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.vcf.gz",
        complete_assembly_phase_assembly_hap2_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.vcf.gz"
    output:
        complete_assembly_phase_assembly_hap1_polish_region_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.polish_region.vcf.gz",
        complete_assembly_phase_assembly_hap2_polish_region_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.polish_region.vcf.gz"
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        bcftools view -R {input.complete_assembly_hap1_polish_region_bed} {input.complete_assembly_phase_assembly_hap1_vcf} -o {output.complete_assembly_phase_assembly_hap1_polish_region_vcf}
        
        bcftools view -R {input.complete_assembly_hap2_polish_region_bed} {input.complete_assembly_phase_assembly_hap2_vcf} -o {output.complete_assembly_phase_assembly_hap2_polish_region_vcf}
        
        tabix {output.complete_assembly_phase_assembly_hap1_polish_region_vcf}
        tabix {output.complete_assembly_phase_assembly_hap2_polish_region_vcf}
        """
        
rule complete_assembly_polish_region_consensus:
    input:
        complete_assembly_phase_assembly_hap1_polish_region_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.phase_assembly.polish_region.vcf.gz",
        complete_assembly_phase_assembly_hap2_polish_region_vcf = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.phase_assembly.polish_region.vcf.gz",
        complete_assembly_hap1_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.fasta"
    output:
        complete_assembly_hap1_polish_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.fasta",
        complete_assembly_hap1_polish_fai = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.fasta.fai",
        complete_assembly_hap2_polish_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.fasta",
        complete_assembly_hap2_polish_fai = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.fasta.fai" 
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        bcftools consensus -H 2 -f {input.complete_assembly_hap1_fa} {input.complete_assembly_phase_assembly_hap1_polish_region_vcf} > {output.complete_assembly_hap1_polish_fa}
        samtools faidx {output.complete_assembly_hap1_polish_fa}
        bcftools consensus -H 2 -f {input.complete_assembly_hap2_fa} {input.complete_assembly_phase_assembly_hap2_polish_region_vcf} > {output.complete_assembly_hap2_polish_fa}
        samtools faidx {output.complete_assembly_hap2_polish_fa}
        """

rule complete_assembly_polish_region_merqury_clip:
    input:
        complete_assembly_hap1_polish_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.fasta",
        complete_assembly_hap1_polish_fai = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.fasta.fai",
        complete_assembly_hap2_polish_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.fasta",
        complete_assembly_hap2_polish_fai = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.fasta.fai",
        sample_WGS_meryl = "c3_merge_snv/meryl/{sample}/{sample}-WGS.meryl/merylIndex"
    output:
        complete_assembly_hap1_polish_clip_region = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.clip.region",
        complete_assembly_hap2_polish_clip_region = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.clip.region",
        final_hap1_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.clip.fasta",
        final_hap2_fa = "c8_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.clip.fasta"
    params:
        sample_WGS_meryl_dir = lambda wildcards, input: os.path.dirname(input.sample_WGS_meryl)
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        mkdir c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.merqury
        cd c8_diploid_path_infer/result/{wildcards.sample}/{wildcards.sample}.merqury
        merqury.sh {params.sample_WGS_meryl_dir} {input.complete_assembly_hap1_polish_fa} {input.complete_assembly_hap2_polish_fa} merqury
        
        bedtools merge -i {wildcards.sample}.hap1.complete_assembly.polish_only.bed -d 2000 | awk '{{if($3-$2>10000)print$0}}' | bedtools complement -i - -g {input.complete_assembly_hap1_polish_fai} | awk '{{if($3-$2>=10000) print $1":"$2+1"-"$3}}' > {output.complete_assembly_hap1_polish_clip_region}
        bedtools merge -i {wildcards.sample}.hap2.complete_assembly.polish_only.bed -d 2000 | awk '{{if($3-$2>10000)print$0}}' | bedtools complement -i - -g {input.complete_assembly_hap2_polish_fai} | awk '{{if($3-$2>=10000) print $1":"$2+1"-"$3}}' > {output.complete_assembly_hap2_polish_clip_region}
        
        samtools faidx {input.complete_assembly_hap1_polish_fa} -r {output.complete_assembly_hap1_polish_clip_region} | awk -v sample={wildcards.sample} 'BEGIN{{sum=1}}{{if(substr($0,1,1)==">"){{split($1,a,">");split(a[2],b,":");print ">"sample"_hap1_ctg"sum"_"b[1];sum+=1}} else print $0}}' > {output.final_hap1_fa}
        
        samtools faidx {input.complete_assembly_hap2_polish_fa} -r {output.complete_assembly_hap2_polish_clip_region} | awk -v sample={wildcards.sample} 'BEGIN{{sum=1}}{{if(substr($0,1,1)==">"){{split($1,a,">");split(a[2],b,":");print ">"sample"_hap2_ctg"sum"_"b[1];sum+=1}} else print $0}}' > {output.final_hap2_fa}
        """
        
# rule xxx:
#     input:
#     output:
#     threads:
#     shell:
#         """
        
#         """
