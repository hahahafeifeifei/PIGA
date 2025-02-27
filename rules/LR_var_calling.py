rule all_LR_var_calling:
    input:
        "c2_call_lr_snv/lr_beagle/concat/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.beagle.filter.vcf.gz"
,


def get_hifi_input_fastqs(wildcards):
    return config["lr_hifi_fastqs"][wildcards.sample]
def get_zmw_input_fastqs(wildcards):
    return config["lr_zmw_fastqs"][wildcards.sample]


rule lr_hifi_pbmm2_map:
    input:
        fastqs = get_hifi_input_fastqs,
        ref = config['reference']['CHM13']
    output:
        hifi_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.hifi.pbmm2.bam"
    params:
        pbmm2 = TOOLS['pbmm2']
    threads: 4
    shell:
        """
        pbmm2 align -j {threads} -J 1 \
            {input.ref} \
            {input.fastqs} \
            {output.hifi_bam} \
            --preset CCS \
            --log-level INFO --sort \
            --rg "@RG\\tID:{wildcards.sample}.Q20\\tSM:{wildcards.sample}.Q20" \
            --sample {wildcards.sample}.Q20
        """

rule lr_zmw_pbmm2_map:
    input:
        fastqs = get_zmw_input_fastqs,
        ref = config['reference']['CHM13']
    output:
        zmw_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.zmw.pbmm2.bam"
    params:
        pbmm2 = TOOLS['pbmm2']
    threads: 4
    shell:
        """
        pbmm2 align -j {threads} -J 1 \
            {input.ref} \
            {input.fastqs} \
            {output.zmw_bam} \
            --log-level INFO --sort \
            --rg "@RG\\tID:{wildcards.sample}.zmw\\tSM:{wildcards.sample}.zmw" \
            --sample {wildcards.sample}.zmw
        """

rule lr_hifi_dv:
    input:
        hifi_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.hifi.pbmm2.bam",
        ref = config['reference']['CHM13']
    output:
        dv_vcf = "c2_call_lr_snv/lr_dv/{sample}/{sample}.deepvariant.vcf.gz",
        dv_gvcf = "c2_call_lr_snv/lr_dv/{sample}/{sample}.deepvariant.gvcf.gz"
    params:
        singularity = TOOLS['singularity'],
        deepvariant = TOOLS['deepvariant']
    shell:
        """
        mkdir -p c2_call_lr_snv/lr_dv/{wildcards.sample}/dv_intermediate_outputs/
        mkdir -p c2_call_lr_snv/lr_dv/{wildcards.sample}/tmp

        REF_DIR={os.path.dirname(input.ref)}
        
        singularity exec -B $(pwd):/project \
            -B $(pwd)/c2_call_lr_snv/lr_dv/{wildcards.sample}/tmp:/tmp \
            -B ${{REF_DIR}}:${{REF_DIR}} \
                ~/software/deepvariant/deepvariant_v1.3_sandbox \
                /opt/deepvariant/bin/run_deepvariant \
                --num_shards {threads} \
                --model_type=PACBIO \
                --ref={input.ref} \
                --reads=/project/{input.hifi_bam} \
                --output_gvcf=/project/{input.dv_gvcf} \
                --output_vcf=/project/{input.dv_vcf} \
                --make_examples_extra_args="vsc_min_count_snps=1,vsc_min_fraction_snps=0.12,vsc_min_count_indels=2,vsc_min_fraction_indels=0.06" \
                --intermediate_results_dir=/project/c2_call_lr_snv/lr_gvcf/{wildcards.sample}/dv_intermediate_outputs/
                --sample_name {wildcards.sample}

        rm -rf c2_call_lr_snv/lr_dv/{wildcards.sample}/dv_intermediate_outputs/*
        """

rule prepare_chr_bed:
    input:
        ref_fai = config['reference']['CHM13'] + '.fai'
    output:
        chr_bed = "c2_call_lr_snv/lr_glnexus/chr_bed/CHM13.{chr}.bed"
    shell:
        """
        awk -v OFS='\\t' -v chr={wildcards.chr} '{{if($1==chr)print $1,"0",$2}}' {input.ref_fai} > {output.chr_bed}
            
        """


rule lr_glnexus:
    input:
        chr_bed = "c2_call_lr_snv/lr_glnexus/chr_bed/CHM13.{chr}.bed",
        dv_gvcfs = expand("c2_call_lr_snv/lr_dv/{sample}/{sample}.deepvariant.gvcf.gz", sample=config['samples'])
    output:
        dv_gvcf_list = "c2_call_lr_snv/lr_glnexus/CKCG.deepvariant.gvcf.list",
        glnexus_bcf = "c2_call_lr_snv/lr_glnexus/chr_result/CKCG.deepvariant.analysis_set.{chr}.bcf",
        glnexus_vcf = "c2_call_lr_snv/lr_glnexus/chr_result/CKCG.deepvariant.analysis_set.{chr}.vcf.gz",
        glnexus_snp_vcf = "c2_call_lr_snv/lr_glnexus/chr_result/CKCG.deepvariant.analysis_set.biallelic_snp.{chr}.vcf.gz"
    threads: 16
    resources:
        #200G
        mem_mb = 200000
    params:
        mem = 200,
        glnexus = TOOLS['glnexus'],
        bcftools = TOOLS['bcftools']
    shell:
        """
        if [ ! -s {output.dv_gvcf_list} ];then \
            ls c2_call_lr_snv/lr_dv/*/*.deepvariant.gvcf.gz > {output.dv_gvcf_list}
        fi
        
        mkdir -p c2_call_lr_snv/lr_glnexus/chr_db/{wildcards.chr}_db
        
        glnexus_cli \
        --config DeepVariant \
        --dir c2_call_lr_snv/lr_glnexus/chr_db/{wildcards.chr}_db \
        --bed c2_call_lr_snv/lr_glnexus/chr_bed/CHM13.{wildcards.chr}.bed \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        $(cat CKCG.deepvariant.gvcf.list) > {output.glnexus_bcf}

        bcftools view --threads {threads} {output.glnexus_bcf} | \
            bcftools norm --threads {threads} -m -any |\
            bcftools view --threads {threads} -i "AQ>40" -o {output.glnexus_vcf}

        tabix -f {output.glnexus_vcf}
        
            bcftools norm --threads {threads} -m +any {output.glnexus_vcf} |\
            bcftools view --threads {threads} -m2 -M2 -v snps -o {output.glnexus_snp_vcf}

        tabix -f {output.glnexus_snp_vcf}
        rm -rf c2_call_lr_snv/lr_glnexus/chr_db/{wildcards.chr}_db
        """

rule glnexus_vcf_merge:
    input:
        glnexus_vcfs = expand("c2_call_lr_snv/lr_glnexus/chr_result/CKCG.deepvariant.analysis_set.{chr}.vcf.gz", chr=chr_list),
        glnexus_snp_vcfs = expand("c2_call_lr_snv/lr_glnexus/chr_result/CKCG.deepvariant.analysis_set.biallelic_snp.{chr}.vcf.gz", chr=chr_list),
        ref = config['reference']['CHM13']

    output:
        merged_glnexus_vcf = "c2_call_lr_snv/lr_glnexus/merge_result/CKCG.deepvariant.analysis_set.vcf.gz",
        merged_glnexus_snp_vcf = "c2_call_lr_snv/lr_glnexus/merge_result/CKCG.deepvariant.analysis_set.biallelic_snp.vcf.gz",
        merged_glnexus_norm_snp_vcf = "c2_call_lr_snv/lr_glnexus/merge_result/CKCG.deepvariant.analysis_set.biallelic_snp.norm.vcf.gz"
    params:
        bcftools = TOOLS['bcftools']
    threads: 28
    shell:
        """
        bcftools concat --threads {threads} -Oz -o {output.merged_glnexus_vcf} $(ls c2_call_lr_snv/lr_glnexus/chr_result/CKCG.deepvariant.analysis_set.chr{{1..22},X,Y,M}.vcf.gz)
        bcftools concat --threads {threads} -Oz -o {output.merged_glnexus_snp_vcf} $(ls c2_call_lr_snv/lr_glnexus/chr_result/CKCG.deepvariant.analysis_set.biallelic_snp.chr{{1..22},X,Y,M}.vcf.gz)

        bcftools norm --threads {threads} -f {input.ref} -o {output.merged_glnexus_norm_snp_vcf}
        """

#run whatshap to genotype the vcf(give likelihood.)
rule lr_whatshap_genotype:
    input:
        merged_glnexus_norm_snp_vcf = "c2_call_lr_snv/lr_glnexus/merge_result/CKCG.deepvariant.analysis_set.biallelic_snp.norm.vcf.gz",
        zmw_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.zmw.pbmm2.bam",
        ref = config['reference']['CHM13']
    output:
        sample_snp_vcf = temp("c2_call_lr_snv/lr_glnexus/{sample}/{sample}.deepvariant.biallelic_snp.vcf.gz"),
        sample_whatshap_snp_vcf = temp("c2_call_lr_snv/lr_whatshap/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.vcf.gz"),
        sample_whatshap_reform_snp_vcf = "c2_call_lr_snv/lr_whatshap/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.reform.vcf.gz"
    params:
        sex = get_sex,
        bcftools = TOOLS['bcftools']
    threads: 16
    shell:
        """
        bcftools view --threads {threads} -I -s {wildcards.sample} {input.merged_glnexus_norm_snp_vcf} -o {output.sample_snp_vcf}

        whatshap genotype -o {output.sample_whatshap_snp_vcf} -r {input.ref} --ignore-read-groups --sample {wildcards.sample} {input.sample_snp_vcf} {input.zmw_bam}
        tabix -f {output.sample_whatshap_snp_vcf}

        bash scripts/snv_phase/whatshap_vcf_reform.sh {output.sample_whatshap_snp_vcf} {params.sex} | bgzip -c > {output.sample_whatshap_reform_snp_vcf}
        tabix -f {output.sample_whatshap_reform_snp_vcf}
        """

rule lr_whatshap_vcf_merge:
    input:
        sample_whatshap_reform_snp_vcfs = expand("c2_call_lr_snv/lr_whatshap/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.reform.vcf.gz", sample = config['samples'])
    output:
        sample_whatshap_reform_snp_vcf_list = "c2_call_lr_snv/lr_whatshap/merge/CKCG.deepvariant.whatshap.vcf.list",
        merge_whatshap_snp_vcf = "c2_call_lr_snv/lr_whatshap/merge/CKCG.deepvariant.whatshap.analysis_set.biallelic_snp.vcf.gz"
    threads: 16
    resources:
        mem_mb = 
    shell:
        """
        ls c2_call_lr_snv/lr_whatshap/*/*.deepvariant.whatshap.biallelic_snp.reform.vcf.gz > {output.sample_whatshap_reform_snp_vcf_list}
        
        bcftools merge \
            --threads {threads} \
            -m none \
            -l {output.sample_whatshap_reform_snp_vcf_list} \
            -Oz \
            -o {output.merge_whatshap_snp_vcf}
        
        tabix -f {output.merge_whatshap_snp_vcf}
        """

rule lr_longshot_genotype:
    input:
        sample_snp_vcf = "c2_call_lr_snv/lr_whatshap/{sample}/{sample}.deepvariant.biallelic_snp.vcf.gz",
        zmw_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.zmw.pbmm2.bam",
        ref = config['reference']['CHM13']
    output:
        sample_longshot_snp_vcf = temp("c2_call_lr_snv/lr_longshot/{sample}/{sample}.deepvariant.longshot.biallelic_snp.vcf.gz"),
        sample_longshot_reform_snp_vcf = "c2_call_lr_snv/lr_longshot/{sample}/{sample}.deepvariant.longshot.biallelic_snp.reform.vcf.gz",
    params:
        sex = get_sex
    threads: 4
    resources:
        #60G
        mem_mb = 60000
    shell:
        """
        longshot \
            --output-ref --force_overwrite \
            -b {input.zmw_bam} \
            -f {input.ref} \
            -v {input.sample_snp_vcf} \
            -s {wildcards.sample} \
            -o {output.sample_longshot_snp_vcf}

        bash scripts/snv_phase/longshot_vcf_reform.sh {output.sample_longshot_snp_vcf} {params.sex} | bgzip -@ {threads} -c  > {output.sample_longshot_reform_snp_vcf}
        """

rule lr_longshot_vcf_merge:
    input:
        sample_longshot_reform_snp_vcfs = expand("c2_call_lr_snv/lr_longshot/{sample}/{sample}.deepvariant.longshot.biallelic_snp.reform.vcf.gz", sample = config['samples'])
    output:
        sample_longshot_reform_snp_vcf_list = "c2_call_lr_snv/lr_longshot/merge/CKCG.deepvariant.longshot.vcf.list",
        merge_longshot_snp_vcf = "c2_call_lr_snv/lr_longshot/merge/CKCG.deepvariant.longshot.analysis_set.biallelic_snp.vcf.gz"
    threads: 28
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        ls c2_call_lr_snv/lr_longshot/*/*.deepvariant.longshot.biallelic_snp.reform.vcf.gz > {output.sample_longshot_reform_snp_vcf_list}
        
        bcftools merge \
            --threads {threads} \
            -m any \
            -l {output.sample_longshot_reform_snp_vcf_list} \
            -Oz \
            -o {output.merge_longshot_snp_vcf}
        
        tabix -f {output.merge_longshot_snp_vcf}
        """

rule lr_genotyped_vcf_gdis_filter:
    input:
        merge_whatshap_snp_vcf = "c2_call_lr_snv/lr_whatshap/merge/CKCG.deepvariant.whatshap.analysis_set.biallelic_snp.vcf.gz",
        merge_longshot_snp_vcf = "c2_call_lr_snv/lr_longshot/merge/CKCG.deepvariant.longshot.analysis_set.biallelic_snp.vcf.gz"
    output:
        gdis_list = "c2_call_lr_snv/gdis_filter/CKCG.deepvariant.whatshap.analysis_set.biallelic_snp.gdis.list.gz",
        whatshap_filter_snp_vcf = "c2_call_lr_snv/gdis_filter/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.vcf.gz"
    threads: 24
    resources:
        mem_mb = 30000
    shell:
        """
        bcftools merge --force-samples -m any \
            --threads {threads} \
            {output.merge_whatshap_snp_vcf} \
            {output.merge_longshot_snp_vcf} | \
            grep -v "^#" | \
            awk -v OFS='\\t' '{{if(substr($0,1,1)=="#")print$0;else {{sum_cor=0; for(i=1;i<=1089;i++) {{split($(i+9),a,":");split($(i+1098),b,":"); if(a[1]==b[1] || (a[1]=="0/1" && b[1]=="1/0") || (a[1]=="1/0" && b[1]=="0/1" || (a[1]=="1/1" && b[1]=="1/1") ) ) sum_cor+=1 }} if(sum_cor<=980) print $1,$2,$4,$5,"GDIS"}} }}' > c2_call_lr_snv/gdis_filter/CKCG.deepvariant.whatshap.analysis_set.biallelic_snp.gdis.list

        bgzip -f c2_call_lr_snv/gdis_filter/CKCG.deepvariant.whatshap.analysis_set.biallelic_snp.gdis.list
        tabix -f -s1 -b2 -e2 {output.gdis_list}

        echo '##FILTER=<ID=GDIS,Description="Genotype concordance between whatshap and longshot is less than 10%">' > c2_call_lr_snv/gdis_filter/header
        
        bcftools annotate --threads {threads} \
            -h c2_call_lr_snv/gdis_filter/header \
            -a {output.gdis_list} \
            -c CHROM,POS,REF,ALT,FILTER \
            {input.merge_whatshap_snp_vcf} | \
            bcftools view --threads {threads} -f . | \
            bcftools plugin fill-tags --threads {threads} | \
            bcftools view --threads {threads} -i "HWE>=1e-6 || ExcHet>=0.1" | bcftools view --threads {threads} -i "(CHR!='chrX' && MEDIAN(GQ)>50 && AC!=0) || (CHR=='chrX' && MEDIAN(GQ)>40 && AC!=0)" -o {output.whatshap_filter_snp_vcf}
        """

# rule lr_re_whatshap:
#     input:
#         whatshap_filter_snp_vcf = "c2_call_lr_snv/gdis_filter/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.vcf.gz",
#         zmw_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.zmw.pbmm2.bam",
#         ref = config['reference']['CHM13'],
#         rename_concat_shapeit4_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.maf0.01.rename.vcf.gz"
#     output:
#         sample_whatshap_filter_snp_vcf = temp("c2_call_lr_snv/lr_re_whatshap/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.vcf.gz"),
#         sample_scaffold = temp("c2_call_lr_snv/lr_re_whatshap/{sample}/{sample}.maf0.01.scaffold.vcf.gz"),
#         sample_whatshap_regenotype_vcf = temp("c2_call_lr_snv/lr_re_whatshap/{sample}/{sample}.deepvariant.whatshap.regenotyping.biallelic_snp.vcf.gz"),
#         sample_whatshap_regenotype_reform_vcf = "c2_call_lr_snv/lr_re_whatshap/{sample}/{sample}.deepvariant.whatshap.regenotyping.biallelic_snp.reform.vcf.gz"
#     params:
#         sex = get_sex
#     threads:
#     resources:
#         mem_mb = 
#     shell:
#         """
#         bcftools view --threads {threads} -I -s {wildcards.sample} {input.whatshap_filter_snp_vcf} -o {output.sample_whatshap_filter_snp_vcf}

#         bcftools view --threads {threads} -I -s {wildcards.sample} {input.rename_concat_shapeit4_vcf} -o {output.sample_scaffold}

#         whatshap genotype \
#             -o {output.sample_whatshap_regenotype_vcf} \
#             -r {input.ref} \
#             --ignore-read-groups \
#             --sample {wildcards.sample} \
#             {output.sample_whatshap_filter_snp_vcf} \
#             {input.zmw_bam} \
#             {output.sample_scaffold}

#         bash scripts/snv_phase/whatshap_vcf_reform.sh {output.sample_whatshap_regenotype_vcf} {params.sex} | bgzip -c > {output.sample_whatshap_regenotype_reform_vcf}
#         """

# rule lr_re_whatshap_vcf_merge:
#     input:
#         sample_whatshap_regenotype_reform_vcfs = expand("c2_call_lr_snv/lr_re_whatshap/{sample}/{sample}.deepvariant.whatshap.regenotyping.biallelic_snp.reform.vcf.gz", sample=config['samples']),
        
#     output:
#         sample_whatshap_regenotype_reform_vcf_list = "c2_call_lr_snv/lr_re_whatshap/merge/vcf.list",
#         merge_whatshap_regenotype_vcf = "c2_call_lr_snv/lr_re_whatshap/merge/CKCG.deepvariant.whatshap.regenotyping.biallelic_snp.vcf.gz",
#         merge_whatshap_regenotype_filter_vcf = "c2_call_lr_snv/lr_re_whatshap/merge/CKCG.deepvariant.whatshap.regenotyping.filter.biallelic_snp.vcf.gz"
#     threads: 16
#     resources:
#         mem_mb = 
#     shell:
#         """
#         ls c2_call_lr_snv/lr_re_whatshap/{wildcards.sample}/*/*.deepvariant.whatshap.regenotyping.biallelic_snp.reform.vcf.gz > {output.sample_whatshap_regenotype_reform_vcf_list}

#         bcftools merge --threads {threads} \
#             -m none \
#             -l {output.sample_whatshap_regenotype_reform_vcf_list} \
#             -Oz \
#             -o {output.merge_whatshap_regenotype_vcf}
        
#         bcftools plugin fill-tags \
#             --threads {threads} \
#             {output.merge_whatshap_regenotype_vcf} | \
#             bcftools view --threads {threads} -i "HWE>=1e-6 || ExcHet>=0.1" | \
#             bcftools view --threads {threads} \
#             -i "(CHR!='chrX' && MEDIAN(GQ)>50 && AC!=0) || (CHR=='chrX' && MEDIAN(GQ)>40 && AC!=0)" \
#             -o {output.merge_whatshap_regenotype_filter_vcf}

#         tabix -f {output.merge_whatshap_regenotype_filter_vcf}
#         """

rule lr_beagle_split_vcf:
    input:
        whatshap_filter_snp_vcf = "c2_call_lr_snv/gdis_filter/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.vcf.gz"
    output:
        chr_num_vcf = "c2_call_lr_snv/lr_beagle/{chr}/CKCG.deepvariant.whatshap.genotyping.filter.biallelic_snp.{chr}.1.vcf.gz"
    threads: 1
    resources:
        #10G
        mem_mb = 10000
    shell:
        """
        mkdir -p c2_call_lr_snv/lr_beagle/{wildcards.chr}

        bcftools view {input.whatshap_filter_snp_vcf} {wildcards.chr} | \
            java -jar ~/software/beagle/beagle4/splitvcf.jar {wildcards.chr} 50000 3000 c2_call_lr_snv/lr_beagle/{wildcards.chr}/CKCG.deepvariant.whatshap.genotyping.filter.biallelic_snp.{wildcards.chr}

        
        """



# checkpoint for counting all of the "chr"-"num" pairs.
checkpoint generate_beagle_split_num_list:
    input:
        chr_num_vcfs = expand("c2_call_lr_snv/lr_beagle/{chr}/CKCG.deepvariant.whatshap.genotyping.filter.biallelic_snp.{chr}.1.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        beagle_split_num_list = "c2_call_lr_snv/lr_beagle/c2_call_lr_snv/lr_beagle/beagle_split_num.list"
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        """
        ls c2_call_lr_snv/lr_beagle/*/CKCG.deepvariant.whatshap.genotyping.filter.biallelic_snp.*.vcf.gz | cut -d "." --output-delimiter=" " -f 7,8  > {output.beagle_split_num_list}
        """


rule lr_beagle:
    input:
        chr_num_vcf = "c2_call_lr_snv/lr_beagle/{chr}/CKCG.deepvariant.whatshap.genotyping.filter.biallelic_snp.{chr}.{num}.vcf.gz",
    output:
        chr_num_beagle_vcf = "c2_call_lr_snv/lr_beagle/{chr}/CKCG.deepvariant.whatshap.genotyping.filter.beagle.biallelic_snp.{chr}.{num}.vcf.gz"
    threads: 32
    resources:
        #80G
        mem_mb = 80000
    shell:
        """
        java -Xmx60g -jar ~/software/beagle/beagle4/beagle.27Jan18.7e1.jar \
            nthreads={threads} \
            gl={input.chr_num_vcf} \
            out=c2_call_lr_snv/lr_beagle/{chr}/CKCG.deepvariant.whatshap.genotyping.filter.beagle.biallelic_snp.{wildcards.chr}.{wildcards.num}

        rm c2_call_lr_snv/lr_beagle/{wildcards.chr}/CKCG.deepvariant.whatshap.genotyping.filter.biallelic_snp.{wildcards.chr}.{wildcards.num}.vcf.gz
        """


def get_chr_beagle_result_files(wildcards, prefix):
    chr_num_combination = checkpoints.generate_beagle_split_num_list.get().output[0]
    pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, num = line.strip().split(" ")
            if chr == wildcards.chr:
                pairs.append(f"c2_call_lr_snv/lr_beagle/{chr}/{prefix}.deepvariant.whatshap.genotyping.filter.beagle.biallelic_snp.{chr}.{num}.vcf.gz")
    return pairs


rule chr_num_beagle_vcf_merge:
    input:
        chr_num_beagle_vcfs = partial(get_chr_beagle_result_files, prefix = config['prefix'])
    output:
        chr_beagle_vcf = "c2_call_lr_snv/lr_beagle/{chr}/CKCG.deepvariant.whatshap.genotyping.filter.beagle.biallelic_snp.{chr}.vcf.gz"
    threads: 1
    resources:
        #20G
        mem_mb = 20000
    shell:
        """
        java -Xmx20g -jar ~/software/beagle/beagle4/mergevcf.jar \
            {wildcards.chr} \
            c2_call_lr_snv/lr_beagle/{wildcards.chr}/CKCG.deepvariant.whatshap.genotyping.filter.beagle.biallelic_snp.{wildcards.chr}.*.vcf.gz | bgzip -c > {output.chr_beagle_vcf}
        tabix -f {output.chr_beagle_vcf}
        """

rule concat_beagle_vcf:
    input:
        chr_beagle_vcf = expand("c2_call_lr_snv/lr_beagle/{chr}/CKCG.deepvariant.whatshap.genotyping.filter.beagle.biallelic_snp.{chr}.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        chr_beagle_vcf_list = "c2_call_lr_snv/lr_beagle/concat/vcf.list",
        concat_beagle_vcf = "c2_call_lr_snv/lr_beagle/concat/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.beagle.vcf.gz",
        concat_beagle_filter_vcf = "c2_call_lr_snv/lr_beagle/concat/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.beagle.filter.vcf.gz"
    threads: 24
    resources:
        #80G
        mem_mb = 80000
    shell:
        """
        
        ls c2_call_lr_snv/lr_beagle/*/CKCG.deepvariant.whatshap.genotyping.filter.beagle.biallelic_snp.*.vcf.gz > {output.chr_beagle_vcf_list}
        
        bcftools concat --threads {threads} \
            -a -f {output.chr_beagle_vcf_list} \
            -Oz \
            -o {output.concat_beagle_vcf}
            
        bcftools plugin fill-tags --threads {threads} {output.concat_beagle_vcf} | bcftools view --threads {threads} -i "AC!=0" -o {output.concat_beagle_filter_vcf}
        """
# rule xxx:
#     input:
#     output:
#     threads:
#     resources:
#         mem_mb = 
#     shell:
#         """
#         """