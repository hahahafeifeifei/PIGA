
chr_list= [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

rule all_infer_diploid_path:
    input:
        expand("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.polish.clip.fasta", sample=config['samples']),
        expand("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.polish.clip.fasta", sample=config['samples'])

rule kmc_prepare_sample_kmer:
    input:
        sr_fq1 = config['sr_fastqs'][0],
        sr_fq2 = config['sr_fastqs'][1],
        hifi_fq = config['lr_hifi_fastqs']
    output:
        kmer_list = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.kmer.list",
        kff = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.kmer.kff"
    params:
        tmp_dir = "c8_diploid_path_infer/sample_assembly/{sample}/tmp"
        prefix = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.kmer"
    threads: 4
    resources:
        mem_mb= 30*1024,
    shell:
        """
        > {output.kmer_list}
        echo {input.sr_fq1} >> {output.kmer_list}
        echo {input.sr_fq2} >> {output.kmer_list}
        echo {input.hifi_fq} >> {output.kmer_list}
        
        mkdir -p {params.tmp_dir}
        kmc -k29 -m28 -o kff -t {threads} -hp @{output.kmer_list} {output.kff} {params.prefix}
        rm -rf {params.tmp_dir}
        """
        
rule get_sample_haplotype_path:
    input:
        hapl = get_hapl_input,
        gbz = get_gbz_input,
        sample_kff = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.kmer.kff"
    output:
        sample_haplotype_gbz = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.haplotypes.gbz"),
        sample_haplotype_gfa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.haplotypes.gfa")
    threads: 4
    resources:
        mem_mb= 450*1024,
    shell:
        """
        vg haplotypes -t {threads} -i {input.hapl} -k {input.sample_kff} -g {output.sample_haplotype_gbz} {input.gbz} --num-haplotypes 17 --include-reference --set-reference {wildcards.sample}
        vg convert -f {output.sample_haplotype_gbz} | \
        awk -v FS='\\t' -v OFS='\\t' '{{if($2=="recombination") {{if($3==1) $2="recombination_ref"; else $3=$3-1}} print $0}}' > {output.sample_haplotype_gfa}
        """
        
rule form_sample_merge_chr_gfa:
    input:
        gfa = get_gfa_input,
        variant_path = get_variant_path_input,
        sample_haplotype_gfa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.haplotypes.gfa"
    output:
        sample_merge_gfa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.merge.gfa"),
        rmac0_gfa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.merge.rmac0.gfa"),
        ref_fa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.ref.fasta")
    wildcard_constraints:
        chr = "|".join(chr_list)
    resources:
        mem_mb= 30*1024,
    shell:
        """
        > {output.sample_merge_gfa}
        cat {input.gfa} >> {output.sample_merge_gfa}
        grep "{wildcards.sample}" {input.variant_path} >> {output.sample_merge_gfa}
        grep -P "{wildcards.chr}\\t" {input.sample_haplotype_gfa} >> {output.sample_merge_gfa}
        python3 scripts/graph-simplification/gfa_remove_ac0_reverse.py {input.sample_merge_gfa} {output.rmac0_gfa}
        python3 scripts/graph-genotyping/gfa_fa.py {output.rmac0_gfa} recombination_ref 1 {output.ref_fa}
        """

rule sample_chr_gfa_deconstruct:
    input:
        rmac0_gfa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.merge.rmac0.gfa"
    output:
        unchop_gfa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.merge.unchop.gfa"),
        unchop_gfaffix_gfa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.merge.unchop.gfaffix.gfa"),
        vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.vcf"),
        norm_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.norm.vcf.gz")
    params:
        sex = get_sex
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 100 * 1024 * attempt,
    shell:
        """
        vg mod -u {input.rmac0_gfa} | vg view - > {output.unchop_gfa}
        gfaffix {output.unchop_gfa} > {output.unchop_gfaffix_gfa} 
        vg deconstruct -t {threads} -P recombination_ref -a -e -C {output.unchop_gfaffix_gfa} > {output.vcf}
        python3 ~/software/script/graph-genotyping/pangenie_norm.py {output.vcf} {wildcards.sample},{wildcards.sample}.snv,recombination {wildcards.sample} {params.sex} | \
        bcftools view -o {output.norm_vcf}
        tabix -f {output.norm_vcf}
        """

rule ref_fasta_vcf_merge:
    input:
        ref_fastas = expand("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.ref.fasta", chr = chr_list),
        norm_chr_vcfs = expand("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.{chr}.norm.vcf.gz", chr = chr_list)
    output:
        merge_ref_fasta = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.fasta",
        merge_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.norm.vcf.gz"),
        vcfbub_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.norm.vcfbub.vcf")
    params:
        sex = get_sex
    shell:
        """
        if [ {params.sex} -eq "female" ];
        then
                cat c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.chr{{{{1..22}},X,M}}.ref.fasta > {output.merge_ref_fasta}
        else
                cat c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.chr{{{{1..22}},X,Y,M}}.ref.fasta > {output.merge_ref_fasta}
        fi
        bcftools concat --threads {threads} {input.norm_chr_vcfs} -o {output.merge_vcf}
        tabix -f {output.merge_vcf}
        vcfbub --input {input.merge_vcf} -l 0 -a 50000 | bcftools view -a - | bcftools view -e "ALT=='.'" -o {output.vcfbub_vcf}
        """


rule vcfbub_vcf_genotype:
    input:
        kmer_list = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.kmer.list",
        vcfbub_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.norm.vcfbub.vcf",
        merge_ref_fasta = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.fasta"
    output:
        kmer_fa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.kmer.fa"),
        genotype_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_genotyping.vcf.gz"),
        phasing_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_phasing.vcf.gz"),
        vcfbub_vcf_gz = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.norm.vcfbub.vcf.gz")
    params:
        pre_prefix = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_pre"
        prefix = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan"
    container:
    resources:
       mem_mb= 120*1024,
    threads: 4
    shell:
        """
        cat {input.kmer_list} | while read line;
        do
            seqkit fq2fa -j {threads} $line
        done > {output.kmer_fa}

        PanGenie-index -v {input.vcfbub_vcf} -r {input.merge_ref_fasta} -t {threads} -o {params.pre_prefix}
        PanGenie -j {threads} -t {threads} -i {input.kmer_fa} -f {params.pre_prefix} -o {params.prefix}
        PanGenie -j {threads} -t {threads} -i {input.kmer_fa} -f {params.pre_prefix} -p -o {params.prefix}
        
        bgzip -f -@ {threads} c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.pan_genotyping.vcf
        tabix -f {output.genotype_vcf}
        
        bgzip -f -@ {threads} c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.pan_phasing.vcf
        tabix -f {output.phasing_vcf}

        bgzip -f -@ {threads} {input.vcfbub_vcf}
        tabix {output.vcfbub_vcf_gz}
        
        rm {params.pre_prefix}* *histo
        """

rule pangenie_refine:
    input:
        genotype_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_genotyping.vcf.gz",
        phasing_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_phasing.vcf.gz",
        vcfbub_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.norm.vcfbub.vcf.gz"
    output:
        phasing_refine_vcf_gz = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_phasing.refine.vcf.gz"),
        genotype_refine_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_genotyping.refine.vcf"),
        genotype_final_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.genotype.vcf.gz"
    params:
        phasing_refine_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_phasing.refine.vcf")
    threads: 4
    resources:
        mem_mb= 30*1024,
    shell:
        """
        python3 scripts/graph-simplification/pangenie_phase_fix.py {input.phasing_vcf} {params.phasing_refine_vcf}
        bgzip -f -@ {threads} {params.phasing_refine_vcf}
        tabix -f {output.phasing_refine_vcf_gz}

        bcftools merge -m all --force-samples {input.genotype_vcf} {input.vcfbub_vcf} | \
            python3 scripts/graph-simplification/graph_genotype_refine.py - {output.genotype_refine_vcf} sample {wildcards.sample}.snv {wildcards.sample}
        bcftools view -a {output.genotype_refine_vcf} | bcftools view -i "ALT!='.'" -o {output.genotype_final_vcf}
        tabix -f {output.genotype_final_vcf}
        """  
        

rule genotype_phase:
    input:
        genotype_final_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.genotype.vcf.gz",
        phasing_refine_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pan_phasing.refine.vcf.gz",
        vcfbub_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.norm.vcfbub.vcf.gz"
    output:
        pangenie_phase_vcf_gz = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pangenie.phase.vcf.gz",
        assembly_phase_vcf_gz = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.assembly.phase.vcf.gz",
        snv_phase_vcf_gz = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.snv.phase.vcf.gz"
    params:
        pangenie_phase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pangenie.phase.vcf",
        assembly_phase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.assembly.phase.vcf",
        snv_phase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.snv.phase.vcf"
    threads: 4
    resources:
        mem_mb= 30*1024,
    shell:
        """
        bcftools merge -m all --force-samples {input.genotype_final_vcf} {input.phasing_refine_vcf} | \
            python3 scripts/graph-genotyping/genotype_phase.py - {params.pangenie_phase_vcf}
        bgzip -@ {threads} -f {params.pangenie_phase_vcf}
        tabix -f {output.pangenie_phase_vcf_gz}
        
        bcftools merge -m all --force-samples {input.genotype_final_vcf} {input.vcfbub_vcf} | \
            bcftools view -s sample,{wildcards.sample} | \
            python3 scripts/graph-genotyping/genotype_phase.py - {params.assembly_phase_vcf}
        bgzip -@ {threads} -f {params.assembly_phase_vcf}
        tabix -f {output.assembly_phase_vcf_gz}

        
        bcftools merge -m all --force-samples {input.genotype_final_vcf} {input.vcfbub_vcf} | \
            bcftools view -s sample,{wildcards.sample}.snv | \
            python3 scripts/graph-genotyping/genotype_phase.py - {params.snv_phase_vcf}
        bgzip -@ {threads} -f {params.snv_phase_vcf}
        tabix -f {output.snv_phase_vcf_gz}
        """  

rule sample_ref_realign:
    input:
        zmw_fq = config['lr_zmw_fastqs'],
        merge_ref_fasta = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.fasta",
    output:
        zmw_bam = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.zmw.bam")
    threads: 8
    resources:
        mem_mb= 60*1024,
    shell:
        """
        minimap2 -t {threads} -ax map-pb -Y -L --eqx --cs {input.merge_ref_fasta} {input.zmw_fq} | \
            samtools view -@ {threads} -hb - | samtools sort -@ {threads} -o {output.zmw_bam}
        samtools index -@ {threads} {output.zmw_bam}
        """

rule hiphase_phase:
    input:
        merge_ref_fasta = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.fasta",
        zmw_bam = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.zmw.bam",
        genotype_final_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.genotype.vcf.gz"
    output:
        hiphase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hiphase.phase.vcf.gz"
    threads: 8
    resources:
        mem_mb= 60*1024,
    shell:
        """
        hiphase -t {threads} \
            --ignore-read-groups \
            --reference {input.merge_ref_fasta} \
            --vcf {input.genotype_final_vcf} \
            --bam {input.zmw_bam} \
            --output-vcf {output.hiphase_vcf} \
            --global-realignment-cputime 600
        
        tabix -f {output.hiphase_vcf}
        """

rule margin_chr_phase:
    input:
        merge_ref_fasta = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.fasta",
        zmw_bam = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.zmw.bam",
        genotype_final_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.genotype.vcf.gz"
    output:
        chr_margin_vcf_gz = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.margin.{chr}.phased.vcf.gz"
    params:
        prefix = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.margin.{chr}",
    shell:
        """
        margin phase \
            {input.zmw_bam} \
            {input.merge_ref_fasta} \
            {input.genotype_final_vcf} \
            test_data/params/allParams.phase_vcf.ont.sv.json \
            -M -t 1 \
            -o {params.prefix} \
            -r {wildcards.chr}
        
        bgzip -f {params.prefix}.phase.vcf
        tabix -f {output.chr_margin_vcf}
        """  

rule merge_margin_phased_vcf:
    input:
        chr_margin_vcfs = expand("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.margin.{chr}.phased.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        merge_margin_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.margin.phase.vcf.gz"
    threads: 4
    shell:
        """
        bcftools concat --threads {threads} {input.chr_margin_vcfs} -o {output.merge_margin_vcf}
        tabix -f {output.merge_margin_vcf}
        
        rm c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.margin.chr*.chunks.csv c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.margin.chr*.phaseset.bed
        rm c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.margin.chr*.phased.vcf.gz c8_diploid_path_infer/sample_assembly/{wildcards.sample}/{wildcards.sample}.margin.chr*.phased.vcf.gz.tbi
        """  
        
rule integrate_phase:
    input:
        genotype_final_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.genotype.vcf.gz",
        pangenie_phase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.pangenie.phase.vcf.gz",
        assembly_phase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.assembly.phase.vcf.gz",
        snv_phase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.snv.phase.vcf.gz",
        hiphase_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hiphase.phase.vcf.gz",
        concat_margin_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.margin.phase.vcf.gz"
    output:
        integrate_phase_log = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.integrate.phase.log"),
        integrate_phase_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.integrate.phase.vcf.gz"),
        integrate_phase_filter_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.integrate.phase.filter.vcf.gz"
    shell:
        """
        bcftools merge -m all --force-samples \
            {input.genotype_final_vcf} \
            {input.snv_phase_vcf} \
            {input.assembly_phase_vcf} \
            {input.pangenie_phase_vcf} \
            {input.hiphase_vcf} \
            {input.concat_margin_vcf} | \
            python3 scripts/graph-genotyping/graph_phase_intergrate.py - {output.integrate_phase_vcf} sample 2:sample,3:sample,4:sample,5:sample,6:sample > {output.integrate_phase_log}
        
        bcftools view -a {output.integrate_phase_vcf} | bcftools view -i "ALT!='.'" -o {output.integrate_phase_filter_vcf}
        tabix -f {output.integrate_phase_filter_vcf}
        """


rule complete_assembly:
    input:
        merge_ref_fasta = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.ref.fasta",
        integrate_phase_filter_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.integrate.phase.filter.vcf.gz"
    output:
        complete_assembly_hap1_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.fasta"
    params:
        sex = get_sex
    resources:
        max_mem_mb = 30
    shell:
        """
        if [ {params.sex} -eq "female" ]; then
            > {output.complete_assembly_hap1_fa}
            samtools faidx {input.merge_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
            samtools faidx {input.merge_ref_fasta} chrX | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
                
            > {output.complete_assembly_hap2_fa}
            samtools faidx {input.merge_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 2 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
            samtools faidx {input.merge_ref_fasta} chr{{X,M}} | \
                bcftools consensus -H 2 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
                
        else 
            > {output.complete_assembly_hap1_fa}
            samtools faidx {input.merge_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
            samtools faidx {input.merge_ref_fasta} chrY | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap1_fa}
                
            > {output.complete_assembly_hap2_fa}
            samtools faidx {input.merge_ref_fasta} chr{{1..22}} | \
                bcftools consensus -H 2 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
            samtools faidx {input.merge_ref_fasta} chr{{X,M}} | \
                bcftools consensus -H 1 {input.integrate_phase_filter_vcf} >> {output.complete_assembly_hap2_fa}
        fi
        """  

        
rule secphase_correct_bam:
    input:
        phase_assembly_hap1_fa = get_hap1_fa_input,
        phase_assembly_hap2_fa = get_hap2_fa_input,
        complete_assembly_hap1_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.fasta"
    output:
        merge_complete_assembly_fa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.merge.complete_assembly.fasta"),
        complete_assembly_phase_assembly_bam = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.merge.complete_assembly.phase_assembly.bam"),
        secphase_log = temp("c8_diploid_path_infer/sample_assembly/{sample}/secphase_out_dir/secphase.out.log"),
        complete_assembly_correct_bam = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.merge.complete_assembly.phase_assembly.correct.bam"),
        phase_assembly_correct_hap1_fa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.phase_assembly.fasta"),
        phase_assembly_correct_hap2_fa = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.phase_assembly.fasta")
    params:
        secphase_dir = "c8_diploid_path_infer/sample_assembly/{sample}/secphase_out_dir"
    threads: 4
    shell:
        """
        > {output.merge_complete_assembly_fa}
        sed "s/chr/H1_chr/g" {input.complete_assembly_hap1_fa} >> {output.merge_complete_assembly_fa}                                                                       
        sed "s/chr/H2_chr/g" {input.complete_assembly_hap2_fa} >> {output.merge_complete_assembly_fa}

        minimap2 -t {threads} -I 8G -ax asm20 -Y -L --eqx --cs \
            {output.merge_complete_assembly_fa} \
            {input.phase_assembly_hap1_fa} {input.phase_assembly_hap2_fa} | \
            samtools sort -@ {threads} -n | samtools view -@ {threads} -hb \
            > {output.complete_assembly_phase_assembly_bam}
        
        secphase -@ {threads} --hifi -i {output.complete_assembly_phase_assembly_bam} -f {output.merge_complete_assembly_fa} -o {params.secphase_dir}
        correct_bam --threads {threads} -i {output.complete_assembly_phase_assembly_bam} -P {output.secphase_log} -o {output.complete_assembly_correct_bam} --primaryOnly -m 0 -a 0
        
        samtools view -@ {threads} {output.complete_assembly_correct_bam} | grep H1 | awk '{{print ">"$1;print$10}}' | seqkit rmdup -o {output.phase_assembly_correct_hap1_fa}
        samtools view -@ {threads} {output.complete_assembly_correct_bam} | grep H2 | awk '{{print ">"$1;print$10}}' | seqkit rmdup -o {output.phase_assembly_correct_hap2_fa}
        """

        
rule get_phase_assembly_sv_vcf:
    input:
        complete_assembly_hap1_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.fasta",
        phase_assembly_correct_hap1_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.phase_assembly.fasta",
        phase_assembly_correct_hap2_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.phase_assembly.fasta"
    output:
        complete_assembly_phase_assembly_hap1_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.phase_assembly.vcf.gz"),
        complete_assembly_phase_assembly_hap2_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.phase_assembly.vcf.gz")
        complete_assembly_phase_assembly_hap1_sv_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}.hap1.complete_assembly.phase_assembly.sv.vcf"),
        complete_assembly_phase_assembly_hap2_sv_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}.hap2.complete_assembly.phase_assembly.sv.vcf")
    threads: 8
    resources:
        mem_mb= 60*1024,
    shell:
        """
        minimap2 -t {threads} -c --cs -x asm20 -B 2 -E 3,1 -O 6,100 \
            {input.complete_assembly_hap1_fa} {input.phase_assembly_correct_hap1_fa} | \
            sort -k6,6 -k8,8n | paftools.js call -f {input.complete_assembly_hap1_fa} -l 5000 -L 5000 - \
            bcftools view -o {output.complete_assembly_phase_assembly_hap1_vcf}
        tabix -f {output.complete_assembly_phase_assembly_hap1_vcf}
        zcat {output.complete_assembly_phase_assembly_hap1_vcf} | \
            awk '{{if(substr($0,1,1)=="#")print$0;else{{if(length($5)-length($4)>=50 || length($4)-length($5)>=50)print$0}} }}' | \
            truvari anno svinfo > {output.complete_assembly_phase_assembly_hap1_sv_vcf}
        
        minimap2 -t {threads} -c --cs -x asm20 -B 2 -E 3,1 -O 6,100 \
            {input.complete_assembly_hap2_fa} {input.phase_assembly_correct_hap2_fa} | \
            sort -k6,6 -k8,8n | paftools.js call -f {input.complete_assembly_hap2_fa} -l 5000 -L 5000 - \
            bcftools view -o {output.complete_assembly_phase_assembly_hap2_vcf}
        tabix -f {output.complete_assembly_phase_assembly_hap2_vcf}
        zcat {output.complete_assembly_phase_assembly_hap2_vcf} | \
            awk '{{if(substr($0,1,1)=="#")print$0;else{{if(length($5)-length($4)>=50 || length($4)-length($5)>=50)print$0}} }}' | \
            truvari anno svinfo > {output.complete_assembly_phase_assembly_hap2_sv_vcf}
        """

   
rule validate_sv:
    input:
        complete_assembly_hap1_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.fasta",
        complete_assembly_phase_assembly_hap1_sv_vcf = "c8_diploid_path_infer/sample_assembly/{sample}.hap1.complete_assembly.phase_assembly.sv.vcf",
        complete_assembly_phase_assembly_hap2_sv_vcf = "c8_diploid_path_infer/sample_assembly/{sample}.hap2.complete_assembly.phase_assembly.sv.vcf",
        zmw_fq = config['lr_zmw_fastqs']
    output:
        complete_assembly_zmw_reads_hap1_bam = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.zmw.bam"),
        complete_assembly_zmw_reads_hap2_bam = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.zmw.bam"),
        complete_assembly_hap1_cutesv_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.phase_assembly.cutesv.vcf"),
        complete_assembly_hap2_cutesv_vcf = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.phase_assembly.cutesv.vcf"),
        complete_assembly_hap1_polish_region_bed = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.phase_assembly.polish_region.bed"),
        complete_assembly_hap2_polish_region_bed = temp("c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.phase_assembly.polish_region.bed")
    params:
        cutesv_tmp = "c8_diploid_path_infer/sample_assembly/{sample}/cutesv.tmp "
    threads: 8
    resources:
        mem_mb= 60*1024,
    shell:
        """
        minimap2 -t {threads} -ax map-pb -Y -L --eqx --cs \
            {input.complete_assembly_hap1_fa} {input.zmw_fq} | \ 
            samtools view -@ {threads} -hb - | \
            samtools sort -@ {threads} -o {output.complete_assembly_zmw_reads_hap1_bam}
        samtools index -@ {threads} {output.complete_assembly_zmw_reads_hap1_bam}

        mkdir -p {params.cutesv_tmp}
        cuteSV {output.complete_assembly_zmw_reads_hap1_bam} {input.complete_assembly_hap1_fa} {output.complete_assembly_hap1_cutesv_vcf} \
            {params.cutesv_tmp} -Ivcf {input.complete_assembly_phase_assembly_hap1_sv_vcf} \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 200 \
            --diff_ratio_merging_DEL 0.5 \
            -q 10 -L -1 -t {threads}
        rm -r {params.cutesv_tmp}

        bcftools view -i "GT!='RR'" {output.complete_assembly_hap1_cutesv_vcf} | \
            bcftools query -f '%CHROM\\t%POS\\t%END\\t%SVLEN\\n' | \
            awk '{{if($4<=0) svlen=-1*$4;else svlen=$4;if(($2-3*svlen)<=0)print $1"\\t0\\t"$2+3*svlen;else print $1"\\t"$2-3*svlen"\\t"$2+3*svlen}}' | \
            bedtools sort -i | bedtools merge -i - > {output.complete_assembly_hap1_polish_region_bed}
        
        
        minimap2 -t {threads} -ax map-pb -Y -L --eqx --cs \
            {input.complete_assembly_hap2_fa} {input.zmw_fq} | \ 
            samtools view -@ {threads} -hb - | \
            samtools sort -@ {threads} -o {output.complete_assembly_zmw_reads_hap2_bam}
        samtools index -@ {threads} {output.complete_assembly_zmw_reads_hap2_bam}

        mkdir -p {params.cutesv_tmp}
        cuteSV {output.complete_assembly_zmw_reads_hap2_bam} {input.complete_assembly_hap2_fa} {output.complete_assembly_hap2_cutesv_vcf} \
            {params.cutesv_tmp} -Ivcf {input.complete_assembly_phase_assembly_hap2_sv_vcf} \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 200 \
            --diff_ratio_merging_DEL 0.5 \
            -q 10 -L -1 -t {threads}
        rm -r {params.cutesv_tmp}
        
        bcftools view -i "GT!='RR'" {output.complete_assembly_hap2_cutesv_vcf} | \
            bcftools query -f '%CHROM\\t%POS\\t%END\\t%SVLEN\\n' | \
            awk '{{if($4<=0) svlen=-1*$4;else svlen=$4;if(($2-3*svlen)<=0)print $1"\\t0\\t"$2+3*svlen;else print $1"\\t"$2-3*svlen"\\t"$2+3*svlen}}' | \
            bedtools sort -i | bedtools merge -i - > {output.complete_assembly_hap2_polish_region_bed}
        """
        

        
rule complete_assembly_polish_region_consensus:
    input:
        complete_assembly_hap1_polish_region_bed = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.phase_assembly.polish_region.bed",
        complete_assembly_hap2_polish_region_bed = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.phase_assembly.polish_region.bed",
        complete_assembly_phase_assembly_hap1_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.phase_assembly.vcf.gz",
        complete_assembly_phase_assembly_hap2_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.phase_assembly.vcf.gz",
        complete_assembly_hap1_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.fasta",
        complete_assembly_hap2_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.fasta"
    output:
        complete_assembly_phase_assembly_hap1_polish_region_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.phase_assembly.polish_region.vcf.gz",
        complete_assembly_phase_assembly_hap2_polish_region_vcf = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.phase_assembly.polish_region.vcf.gz",
        complete_assembly_hap1_polish_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.polish.fasta",
        complete_assembly_hap1_polish_fai = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.polish.fasta.fai",
        complete_assembly_hap2_polish_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.polish.fasta",
        complete_assembly_hap2_polish_fai = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.polish.fasta.fai"
    threads: 4
    resources:
        mem_mb= 30*1024,
    shell:
        """
        bcftools view -R {input.complete_assembly_hap1_polish_region_bed} {input.complete_assembly_phase_assembly_hap1_vcf} -o {output.complete_assembly_phase_assembly_hap1_polish_region_vcf}
        tabix {output.complete_assembly_phase_assembly_hap1_polish_region_vcf}
        bcftools consensus -H 2 -f {input.complete_assembly_hap1_fa} {input.complete_assembly_phase_assembly_hap1_polish_region_vcf} > {output.complete_assembly_hap1_polish_fa}
        samtools faidx {output.complete_assembly_hap1_polish_fa}

        bcftools view -R {input.complete_assembly_hap2_polish_region_bed} {input.complete_assembly_phase_assembly_hap2_vcf} -o {output.complete_assembly_phase_assembly_hap2_polish_region_vcf}
        tabix {output.complete_assembly_phase_assembly_hap2_polish_region_vcf}
        bcftools consensus -H 2 -f {input.complete_assembly_hap2_fa} {input.complete_assembly_phase_assembly_hap2_polish_region_vcf} > {output.complete_assembly_hap2_polish_fa}
        samtools faidx {output.complete_assembly_hap2_polish_fa}
        """


rule complete_assembly_polish_region_merqury_clip:
    input:
        complete_assembly_hap1_polish_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.polish.fasta",
        complete_assembly_hap1_polish_fai = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.polish.fasta.fai",
        complete_assembly_hap2_polish_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.polish.fasta",
        complete_assembly_hap2_polish_fai = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.polish.fasta.fai",
        sample_meryl = get_sample_meryl_input
    output:
        complete_assembly_hap1_polish_clip_region = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.polish.clip_region.bed",
        complete_assembly_hap2_polish_clip_region = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.polish.clip_region.bed",
        final_hap1_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap1.complete_assembly.polish.clip.fasta",
        final_hap2_fa = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.hap2.complete_assembly.polish.clip.fasta"
    params:
        merqury_dir = "c8_diploid_path_infer/sample_assembly/{sample}/{sample}.merqury"
    threads: 4
    resources:
        mem_mb= 30*1024,
    shell:
        """
        mkdir -p {params.merqury_dir}
        cd {params.merqury_dir}
        merqury.sh {params.merqury_dir} {input.complete_assembly_hap1_polish_fa} {input.complete_assembly_hap2_polish_fa} merqury
        
        cd -
        bedtools merge -i {params.merqury_dir}/{wildcards.sample}.hap1.complete_assembly.polish_only.bed -d 2000 | \
            awk '{{if($3-$2>10000)print$0}}' | \
            bedtools complement -i - -g {input.complete_assembly_hap1_polish_fai} | \
            awk '{{if($3-$2>=10000) print $1":"$2+1"-"$3}}' > {output.complete_assembly_hap1_polish_clip_region}
        
        samtools faidx {input.complete_assembly_hap1_polish_fa} -r {output.complete_assembly_hap1_polish_clip_region} | \
        awk -v sample={wildcards.sample} 'BEGIN{{sum=1}}{{if(substr($0,1,1)==">"){{split($1,a,">");split(a[2],b,":");print ">"sample"_hap1_ctg"sum"_"b[1];sum+=1}} else print $0}}' > {output.final_hap1_fa}

        bedtools merge -i {params.merqury_dir}/{wildcards.sample}.hap2.complete_assembly.polish_only.bed -d 2000 | \
            awk '{{if($3-$2>10000)print$0}}' | \
            bedtools complement -i - -g {input.complete_assembly_hap2_polish_fai} | \
            awk '{{if($3-$2>=10000) print $1":"$2+1"-"$3}}' > {output.complete_assembly_hap2_polish_clip_region}
        
        samtools faidx {input.complete_assembly_hap2_polish_fa} -r {output.complete_assembly_hap2_polish_clip_region} | \
        awk -v sample={wildcards.sample} 'BEGIN{{sum=1}}{{if(substr($0,1,1)==">"){{split($1,a,">");split(a[2],b,":");print ">"sample"_hap2_ctg"sum"_"b[1];sum+=1}} else print $0}}' > {output.final_hap2_fa}
        rm -rf {params.merqury_dir}
        """
        
# rule xxx:
#     input:
#     output:
#     threads:
#     shell:
#         """
        
#         """
