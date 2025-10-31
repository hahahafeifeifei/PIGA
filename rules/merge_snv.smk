rule all_merge_snv:
    input:
        f"c3_merge_snv/merged_vcf/{config['prefix']}.consensus.merfin.vcf.gz"

rule sr_lr_bcftools_isec:
    input:
        sr_callset_vcf = get_sr_vcf_input
        lr_callset_vcf = get_lr_vcf_input
    output:
        sr_spec_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_specific.vcf",
        lr_spec_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_specific.vcf",
        sr_shared_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_shared.vcf",
        lr_shared_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_shared.vcf"
    params:
        prefix = f"c3_merge_snv/merged_vcf/srs_lrs_compare"
    threads: 16
    shell:
        """
        bcftools isec --threads {threads} \
            {input.sr_callset_vcf} \
            {input.lr_callset_vcf} \
            -p {params.prefix}

        cp {params.prefix}/0000.vcf {output.sr_spec_vcf}
        cp {params.prefix}/0001.vcf {output.lr_spec_vcf}
        cp {params.prefix}/0002.vcf {output.sr_shared_vcf}
        cp {params.prefix}/0003.vcf {output.lr_shared_vcf}
        """

rule generate_lrs_sub_lists:
    input:
        sr_shared_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_shared.vcf",
        lr_shared_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_shared.vcf"
    output:
        sr_hwe = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_shared.hwe",
        lr_hwe = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_shared.hwe",
        info = f"c3_merge_snv/merged_vcf/{config['prefix']}.shared.info",
        hwe_info = f"c3_merge_snv/merged_vcf/{config['prefix']}.shared.hwe.info",
        lrs_sub_list = f"c3_merge_snv/merged_vcf/{config['prefix']}.shared.sub.list.gz"
    threads: 16

    shell:
        """
        bcftools query {input.lr_shared_vcf} -f "%INFO/HWE\\n" > {output.lr_hwe}
        bcftools query {input.sr_shared_vcf} -f "%INFO/HWE\\n" > {output.sr_hwe}
        bcftools query {input.sr_shared_vcf} -f "%CHROM\\t%POS\\t%REF\\t%ALT\\n" > {output.info}

        paste {output.info} {output.sr_hwe} {output.lr_hwe} | awk -v OFS='\\t' \
        '{{ \
            if ($5 == 0 && $6 == 0) {{ \
                ratio = 0; \
            }} else if ($5 == 0) {{ \
                ratio = 50; \
            }} else if ($6 == 0) {{ \
                ratio = -50; \
            }} else {{ \
                ratio = -1 * log($5 / $6) / log(10); \
            }} \
            print $1, $2, $3, $4, ratio; \
        }}' | \
        awk -v OFS='\\t' '{{if($5>=3)print $1,$2,$3,$4,"LRS_SUB"}}' | \
        sort -k1,1V -k2,2n | bgzip -@ {threads} -c > {output.hwe_info}
        tabix -s1 -b2 -e2 {output.lrs_sub_list}
        """


# based on the lrs_sub_list, remove those variants in sr_callset and select out those variants in lr_callset.
rule lrs_sub_process:
    input:
        sr_shared_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_shared.vcf",
        lr_shared_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_shared.vcf",
        lrs_sub_list = f"c3_merge_snv/merged_vcf/{config['prefix']}.shared.sub.list.gz"
    output:
        header = "c3_merge_snv/merged_vcf/header",
        sr_filter_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_shared.sub_filter.vcf.gz",
        lr_select_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_shared.sub_select.vcf.gz"
    threads: 16
    shell:
        """
        echo '##FILTER=<ID=LRS_SUB,Description="The short-read genotype of variants is substituted by long-read genotype">' > {output.header}
        
        bcftools annotate --threads {threads} -h {output.header}\
            -a {input.lrs_sub_list} \
            -c CHROM,POS,REF,ALT,FILTER \
            {input.sr_shared_vcf} | \
            bcftools view --threads {threads} -e "FILTER=='LRS_SUB'" | \
            bgzip -@ {threads} -c > {output.sr_filter_vcf}
        
        bcftools annotate --threads {threads} -h {output.header} \
            -a {input.lrs_sub_list} \
            -c CHROM,POS,REF,ALT,FILTER \
            {input.lr_shared_vcf} | \
            bcftools view --threads {threads} -i "FILTER=='LRS_SUB'" | \
            sed 's/LRS_SUB/PASS/g' | bgzip -@ {threads} -c > {output.lr_select_vcf}
        """


rule final_merge:
    input:
        sr_filter_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_shared.sub_filter.vcf.gz",
        lr_select_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_shared.sub_select.vcf.gz",
        sr_spec_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.srs_specific.vcf",
        lr_spec_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.lrs_specific.vcf"
    output:
        consensus_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.consensus.vcf.gz"
    threads: 16
    shell:
        """
        bcftools concat {input.sr_filter_vcf} {input.lr_select_vcf} {input.sr_spec_vcf} {input.lr_spec_vcf} \
            -a --threads {threads} -o {output.consensus_vcf}
        tabix -f {output.consensus_vcf}
        """


# merfin polish
#use merylIndex to represent the whole {sample}.meryl directory.
rule prepare_sample_kmer:
    input:
        sr_fq1 = config['sr_fastqs'][0],
        sr_fq2 = config['sr_fastqs'][1],
        hifi_fq = config['lr_hifi_fastqs'],
    output:
        sample_meryl = directory("c3_merge_snv/sample_meryl/{sample}/{sample}.meryl")
    threads: 4
    resources:
        max_mem_gb = 30
    shell:
        """
        meryl count k=21 memory={params.max_mem_gb} threads={threads} output {output.sample_meryl} {input.sr_fq1} {input.sr_fq2} {input.hifi_fq}
        """

rule prepare_ref_kmer:
    input:
        ref = config['reference']['CHM13']
    output:
        ref_meryl = directory("c3_merge_snv/sample_meryl/chm13/chm13.meryl")
    threads: 4
    resources:
        max_mem_gb = 30
    shell:
        """
        meryl count k=21 memory={params.max_mem_gb} threads={threads} output {output.ref_meryl} {output.ref}
        """


rule merfin_filter:
    input:
        consensus_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.consensus.vcf.gz",
        sample_meryl = "c3_merge_snv/sample_meryl/{sample}/{sample}.meryl",
        ref = config['reference']['CHM13'],
        ref_meryl = "c3_merge_snv/sample_meryl/chm13/chm13.meryl"
    output:
        sample_consensus_vcf = "c3_merge_snv/sample_vcf/{sample}/{sample}.consensus.vcf",
        sample_consensus_filter_vcf = "c3_merge_snv/sample_vcf/{sample}/{sample}.consensus.filter.vcf.gz"
    params:
        prefix = "c3_merge_snv/sample_vcf/{sample}/{sample}.consensus"
    threads: 4
    resources:
        max_mem_gb = 60
    shell:
        """
        bcftools view -s {wildcards.sample} --threads {threads} {input.consensus_vcf} | \
            bcftools view --threads {threads} -e 'GT=="0/0" || GT=="0|0" || GT=="0"' | \
            awk -v OFS='\\t' '{{if(substr($1,1,1)=="#") print$0; else{{if(length($10)==1) {{gt=$10"/"$10; print $1,$2,$3,$4,$5,$6,$7,$8,$9,gt }} else print$0}} }}' | \
            bcftools view --threads {threads} -i "GT=='./.' || GT=='.|.'" -o {output.sample_consensus_vcf}
        
        merfin -filter \
            -sequence {input.ref} \
            -readmers {input.sample_meryl} \
            -seqmers {input.ref_meryl} \
            -vcf {output.sample_consensus_vcf} \
            -output {params.prefix} \
            -threads {threads} \
            -memory {params.max_mem_gb}
        bgzip -f -@ {threads} {params.prefix}.filter.vcf
        tabix -f {output.sample_consensus_filter_vcf}
        """

rule merge_merfin_vcf:
    input:
        filter_vcfs = expand("c3_merge_snv/sample_vcf/{sample}/{sample}.consensus.filter.vcf.gz", sample = samples_list)
    output:
        merge_merfin_vcf = f"c3_merge_snv/merged_vcf/{config['prefix']}.consensus.merfin.vcf.gz"
    threads: 16
    resources:
        max_mem_gb = 60
    shell:
        """
        bcftools merge -0 -m none --threads {threads} {input.filter_vcfs} | \
            bcftools sort -m {params.max_mem_gb}G | \
            bcftools plugin fill-tags --threads {threads} | \
            bcftools view --threads {threads} -i "AC!=0 && HWE>=1e-6" -o {output.merge_merfin_vcf}
        tabix {output.merge_merfin_vcf}
        """
