rule all_LR_var_calling:
    input:
        f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.beagle.vcf.gz"

chr_list = [f'chr{i}' for i in range(1, 23)] + ['chrX']

rule lr_hifi_pbmm2_map:
    input:
        hifi_fq = config['lr_hifi_fastqs'],
        ref = config['reference']['CHM13']
    output:
        hifi_bam = "c2_call_lr_snv/sample_bam/{sample}/{sample}.hifi.srt.bam"
    threads: 4
    shell:
        """
        pbmm2 align -j {threads} -J 1 \
            {input.ref} \
            {input.hifi_fq} \
            {output.hifi_bam} \
            --preset CCS \
            --log-level INFO --sort \
            --rg "@RG\\tID:{wildcards.sample}.hifi\\tSM:{wildcards.sample}.hifi" \
            --sample {wildcards.sample}.hifi
        """

rule lr_zmw_pbmm2_map:
    input:
        zmw_fq = config['lr_zmw_fastqs'],
        ref = config['reference']['CHM13']
    output:
        zmw_bam = "c2_call_lr_snv/sample_bam/{sample}/{sample}.zmw.srt.bam"
    threads: 4
    shell:
        """
        pbmm2 align -j {threads} -J 1 \
            {input.ref} \
            {input.zmw_fq} \
            {output.zmw_bam} \
            --log-level INFO --sort \
            --rg "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            --sample {wildcards.sample}
        """

#How to run
rule lr_hifi_dv:
    input:
        hifi_bam = "c2_call_lr_snv/sample_bam/{sample}/{sample}.hifi.srt.bam",
        ref = config['reference']['CHM13']
    output:
        dv_gvcf = "c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.g.vcf.gz"
    threads: 4
    shell:
        """
        mkdir -p c2_call_lr_snv/sample_gvcf/{wildcards.sample}/dv_intermediate_outputs/

        /opt/deepvariant/bin/run_deepvariant \
            --num_shards {threads} \
            --model_type=PACBIO \
            --ref={input.ref} \
            --reads={input.hifi_bam} \
            --output_gvcf={output.dv_gvcf} \
            --make_examples_extra_args="vsc_min_count_snps=1,vsc_min_fraction_snps=0.12,vsc_min_count_indels=2,vsc_min_fraction_indels=0.06" \
            --sample_name {wildcards.sample}

        rm -rf c2_call_lr_snv/sample_gvcf/{wildcards.sample}/dv_intermediate_outputs/*
        """

rule prepare_chr_bed:
    input:
        ref_fai = config['reference']['CHM13'] + '.fai'
    output:
        chr_bed = "c2_call_lr_snv/chr_vcf/chm13.{chr}.bed"
    shell:
        """
        awk -v OFS='\\t' -v chr={wildcards.chr} '{{if($1==chr)print $1,"0",$2}}' {input.ref_fai} > {output.chr_bed}
        """

rule lr_glnexus:
    input:
        ref = config['reference']['CHM13'],
        chr_bed = "c2_call_lr_snv/chr_vcf/chr_bed/chm13.{chr}.bed",
        dv_gvcfs = expand("c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.g.vcf.gz", sample=samples_list)
    output:
        glnexus_bcf = f"c2_call_lr_snv/chr_vcf/{config['prefix']}.deepvariant.{chr}.bcf",
        glnexus_vcf = f"c2_call_lr_snv/chr_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.{chr}.vcf.gz",
        glnexus_vcf_tbi = f"c2_call_lr_snv/chr_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.{chr}.vcf.gz.tbi"
    params:
        chr_db = "c2_call_lr_snv/chr_vcf/{chr}_db"
    threads: 16
    resources:
        max_mem_gb = 200,
    shell:
        """
        ls {input.dv_gvcfs} > {output.gvcf_map}

        mkdir {params.chr_db}
        glnexus_cli \
            --config DeepVariant \
            --dir {params.chr_db} \
            --bed {input.chr_bed} \
            --threads {threads} \
            --mem-gbytes {resources.max_mem_gb} \
            {input.dv_gvcfs} > {output.glnexus_bcf}

        bcftools view --threads {threads} {output.glnexus_bcf} | \
            bcftools norm --threads {threads} -m -any -f {input.ref} |\
            bcftools view --threads {threads} -i "AQ>40" |\
            bcftools norm --threads {threads} -m +any |\
            bcftools view --threads {threads} -m2 -M2 -v snps |\
            bcftools norm --threads {threads} -m -any ｜\
            bcftools sort -o {output.glnexus_vcf}

        tabix -f {output.glnexus_vcf}
        rm -rf {params.chr_db}
        """

rule glnexus_vcf_merge:
    input:
        glnexus_vcfs = expand(f"c2_call_lr_snv/chr_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.{chr}.vcf.gz", chr=chr_list)
    output:
        merged_glnexus_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.vcf.gz"
    threads: 16
    shell:
        """
        bcftools concat --threads {threads} -Oz -o {output.merged_glnexus_vcf} {input.glnexus_vcfs}
        """

rule lr_whatshap_genotype:
    input:
        merged_glnexus_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.vcf.gz",
        zmw_bam = "c2_call_lr_snv/sample_bam/{sample}/{sample}.zmw.srt.bam",
        ref = config['reference']['CHM13']
    output:
        whatshap_snp_vcf = "c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.vcf.gz",
        whatshap_reform_snp_vcf = "c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.reform.vcf.gz"
    params:
        sex = get_sex
    threads: 16
    shell:
        """
        whatshap genotype -o {output.whatshap_snp_vcf} -r {input.ref} --ignore-read-groups --sample {wildcards.sample} {input.merged_glnexus_vcf} {input.zmw_bam}
        tabix -f {output.whatshap_snp_vcf}

        bash scripts/snv_phase/whatshap_vcf_reform.sh {output.whatshap_snp_vcf} {params.sex} | bgzip -c > {output.whatshap_reform_snp_vcf}
        tabix -f {output.whatshap_reform_snp_vcf}
        """

rule lr_whatshap_vcf_merge:
    input:
        whatshap_reform_snp_vcfs = expand("c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.reform.vcf.gz", sample = samples_list)
    output:
        merge_whatshap_snp_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.vcf.gz",
        merge_whatshap_filter_snp_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.vcf.gz",
        merge_whatshap_filter_snp_vcf_tbi = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.vcf.gz.tbi"
    threads: 16
    resources:
        max_mem_gb = 300
    shell:
        """
        bcftools merge --threads {threads} -m none {input.whatshap_reform_snp_vcfs} -Oz -o {output.merge_whatshap_snp_vcf}
        bcftools plugin fill-tags --threads {threads} {output.merge_whatshap_snp_vcf} | \
        bcftools view --threads {threads} -e "(HWE<1e-6 && ExcHet<0.1) || F_MISSING>0.1" | \
        bcftools view --threads {threads} -i "(CHR!='chrX' && MEDIAN(GQ)>50 && AC!=0) || (CHR=='chrX' && MEDIAN(GQ)>40 && AC!=0)" -o {output.merge_whatshap_filter_snp_vcf}
        tabix {output.merge_whatshap_filter_snp_vcf}
        """

rule lr_beagle_split_vcf:
    input:
        whatshap_filter_snp_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.vcf.gz"
    output:
        chr_num_vcf = f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.{chr}.1.vcf.gz"
    params:
        prefix = f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.{chr}"
    threads: 1
    shell:
        """
        mkdir -p c2_call_lr_snv/lr_beagle/{wildcards.chr}

        bcftools view {input.whatshap_filter_snp_vcf} {wildcards.chr} | \
            java -jar splitvcf.jar {wildcards.chr} 50000 3000 {params.prefix}
        
        """

# checkpoint for counting all of the "chr"-"num" pairs.
checkpoint generate_beagle_split_num_list:
    input:
        chr_num_vcfs = expand(f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.{chr}.1.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        beagle_split_num_list = "c2_call_lr_snv/chr_vcf/beagle_split_num.list"
    params:
        prefix = f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter"
    threads: 1
    shell:
        """
        ls {params.prefix}.*.vcf.gz > {output.beagle_split_num_list}
        """


def get_chr_beagle_result_files(wildcards, prefix):
    chr_num_combination = checkpoints.generate_beagle_split_num_list.get().output[0]
    pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, num = line.strip().split("/")[-1].split(".")[5:7]
            if chr == wildcards.chr:
                pairs.append(f"c2_call_lr_snv/chr_vcf/{chr}/{prefix}.deepvariant.whatshap.beagle.{chr}.{num}.vcf.gz")
    return pairs

rule lr_beagle:
    input:
        chr_num_vcf = f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.{chr}.{num}.vcf.gz",
    output:
        chr_num_beagle_vcf = f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.beagle.{chr}.{num}.vcf.gz"
    params:
        prefix = f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.beagle.{chr}.{num}"
    threads: 8
    resources:
        mem_max_gb = 80
    shell:
        """
        java "-Xmx{resources.max_mem_gb}G" -jar beagle.27Jan18.7e1.jar \
            nthreads={threads} \
            gl={input.chr_num_vcf} \
            out={params.prefix}
        """

rule chr_num_beagle_vcf_merge:
    input:
        chr_num_beagle_vcfs = partial(get_chr_beagle_result_files, prefix = config['prefix'])
    output:
        chr_beagle_vcf = f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.beagle.{chr}.vcf.gz"
    threads: 1
    resources:
        mem_max_gb = 20
    shell:
        """
        java "-Xmx{resources.max_mem_gb}G" -jar mergevcf.jar {wildcards.chr} {input.chr_num_beagle_vcfs} | bgzip -c > {output.chr_beagle_vcf}
        tabix -f {output.chr_beagle_vcf}
        """


rule concat_beagle_vcf:
    input:
        chr_beagle_vcfs = expand(f"c2_call_lr_snv/chr_vcf/{chr}/{config['prefix']}.deepvariant.whatshap.beagle.{chr}.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        concat_beagle_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.beagle.vcf.gz"
    resources:
        mem_max_gb = 80
    shell:
        """        
        bcftools concat --threads {threads} -a {input.chr_beagle_vcfs} | \
        bcftools annotate -x FORMAT/DS,FORMAT/GP | \
        bcftools plugin fill-tags --threads {threads} | \
        bcftools view --threads {threads} -i "AC!=0" -o {output.concat_beagle_vcf}
        tabix {output.concat_beagle_vcf}
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
