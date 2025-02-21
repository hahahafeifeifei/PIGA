sr_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.analysis_set.biallelic.vcf.gz"
lr_vcf = "c2_call_lr_snv/lr_beagle/concat/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.beagle.filter.vcf.gz"

rule all:
    input:
        "srs_lrs_compare/CKCG.analysis_set.call_set.consensus.vcf.gz",
        "srs_lrs_compare/CKCG.analysis_set.call_set.consensus.srs_scaffold.info.gz",
        "srs_lrs_compare/CKCG.analysis_set.call_set.consensus.shared.vcf.gz",

# 初始处理步骤

rule process_sr:
    input:
        sr = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.analysis_set.biallelic.vcf.gz"

    output:
        sr_callset_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.srs.vcf.gz"
    threads: 16

    shell:
        """
        bcftools view --threads {threads} {input.sr} | \
            bcftools view --threads {threads} -i "AC!=0" | \
            bcftools annotate --threads {threads} -x QUAL,INFO,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PGT,FORMAT/PID,FORMAT/PL,FORMAT/PS |
            sed 's/-WGS//g' | bgzip -@ {threads} -c > {output.sr_callset_vcf}
        """

rule process_lr:
    input:
        lr = "c2_call_lr_snv/lr_beagle/concat/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.beagle.filter.vcf.gz"

    output:
        lr_callset_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.lrs.vcf.gz"
    threads: 16

    shell:
        """
        zcat {input.lr} | sed 's/-CLR//g' | sed 's/-mix//g' | sed 's/-ccs//g' | bgzip -@ {threads} -c > {output.lr_callset_vcf}
        tabix -f {output.lr_callset_vcf}
        """

rule split_srs:
    input:
        sr_callset_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.srs.vcf.gz"
    output:
        sr_callset_snp_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.srs.snp.vcf.gz",
        sr_callset_indel_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.srs.indel.vcf.gz"
    threads: 16

    shell:
        """
        bcftools view --threads {threads} -v snps {input.sr_callset_vcf} -o {output.sr_callset_snp_vcf}
        tabix -f {output.sr_callset_snp_vcf}
        bcftools view --threads {threads} -V snps {input.sr_callset_vcf} -o {output.sr_callset_indel_vcf}
        tabix -f {output.sr_callset_indel_vcf}
        """


rule sr_lr_bcftools_isec:
    input:
        sr_callset_snp_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.srs.snp.vcf.gz",
        lr_callset_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.lrs.vcf.gz"
    output:
        sr_spec_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_specific.vcf",
        lr_spec_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_specific.vcf",
        sr_shared_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.vcf",
        lr_shared_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.vcf"
    threads: 40

    shell:
        """
        bcftools isec --threads {threads} \
            {input.srs_snp} \
            {input.lrs} \
            -p c3_merge_snv/callset/srs_lrs_compare

        mv c3_merge_snv/callset/srs_lrs_compare/0000.vcf {output.sr_spec_vcf}
        mv c3_merge_snv/callset/srs_lrs_compare/0001.vcf {output.lr_spec_vcf}
        mv c3_merge_snv/callset/srs_lrs_compare/0002.vcf {output.sr_shared_vcf}
        mv c3_merge_snv/callset/srs_lrs_compare/0003.vcf {output.lr_shared_vcf}
        
        """

# HWE分析

rule calc_hwe:
    input:
        sr_shared_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.vcf",
        lr_shared_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.vcf"
    output:
        sr_hwe = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.hwe",
        lr_hwe = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.hwe",
        info = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.info",
        hwe_info = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.hwe.info"
    threads: 16

    shell:
        """
        bcftools query {input.lr_shared_vcf} -f "%INFO/HWE\n" > {output.lr_hwe}
        bcftools plugin fill-tags --threads {threads} {input.sr_shared_vcf} | \
            bcftools query -f "%INFO/HWE\n" > {output.sr_hwe}
        bcftools query {input.sr_shared_vcf} -f "%CHROM\t%POS\t%REF\t%ALT\n" > {output.info}

        paste {input.info} {input.sr_hwe} {input.lr_hwe} | awk -v OFS='\\t' '
            {{if($5==0 && $6==0) ratio=0; 
            else {if($5==0) ratio=50;
            else {if($6==0) ratio=-50; 
            else ratio=-1*log($5/$6)/log(10)}};
            print $1,$2,$3,$4,ratio}}' > {output.hwe_info}
        """


# 生成筛选列表

rule generate_lrs_sub_lists:
    input:
        sr_shared_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.vcf",
        hwe_info = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.hwe.info"
    output:
        lrs_sub_list = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.sub.list.gz"
    threads: 16

    shell:
        """
        awk -v OFS='\\t' '{{if($5>=3)print $1,$2,$3,$4,"LRS_SUB"}}' {input.hwe_info} > c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.hwe_sub.list

        bcftools plugin fill-tags --threads {threads} {input.sr_shared_vcf} | 
        bcftools view -H -i "NS<=1007" | 
        awk -v OFS='\\t' '{{print $1,$2,$4,$5,"LRS_SUB"}}' > c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.miss_sub.list

        cat c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.hwe_sub.list c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.miss_sub.list | sort -k1,1V -k2,2n | bgzip -@ {threads} -c > {output.lrs_sub_list}
        tabix -s1 -b2 -e2 {output.lrs_sub_list}
        """


# based on the lrs_sub_list, remove those variants in sr_callset and select out those variants in lr_callset.
rule lrs_sub_process:
    input:
        sr_shared_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.vcf",
        lr_shared_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.vcf",
        lrs_sub_list = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.shared.sub.list.gz"
    output:
        sr_filter_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.sub_filter.vcf.gz",
        lr_select_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.sub_select.vcf.gz",
    threads: 16

    shell:
        """
        echo '##FILTER=<ID=LRS_SUB,Description="The short-read genotype of variants is substituted by long-read genotype">' > c3_merge_snv/callset/srs_lrs_compare/header
        
        bcftools annotate --threads {threads} -h c3_merge_snv/callset/srs_lrs_compare/header \
            -a {input.lrs_sub_list} \
            -c CHROM,POS,REF,ALT,FILTER \
            {input.sr_shared_vcf} | \
            bcftools view --threads {threads} -e "FILTER=='LRS_SUB'" | \
            bgzip -@ {threads} -c > {output.sr_filter_vcf}
        
        bcftools annotate --threads {threads} -h c3_merge_snv/callset/srs_lrs_compare/header \
            -a {input.lrs_sub_list} \
            -c CHROM,POS,REF,ALT,FILTER \
            {input.lr_shared_vcf} | \
            bcftools view --threads {threads} -i "FILTER=='LRS_SUB'" | \
            sed 's/LRS_SUB/PASS/g' | bgzip -@ {threads} -c > {output.lr_select_vcf}
        """


rule anno_merge:
    input:
        sr_filter_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.sub_filter.vcf.gz",
        lr_select_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.sub_select.vcf.gz",
        sr_spec_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_specific.vcf",
        lr_spec_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_specific.vcf",
        sr_callset_indel_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.srs.indel.vcf.gz"
    output:
        sr_filter_anno = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.sub_filter.anno",
        lr_select_anno = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.sub_select.anno",
        sr_spec_anno = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_specific.anno",
        lr_spec_anno = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_specific.anno",
        sr_indel_anno = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs.indel.anno",
        consensus_anno = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.anno.gz"
    threads: 16
    shell:
        """
        bcftools view --threads {threads} -H {input.sr_filter_vcf} | awk -v OFS='\\t' '{{print$1,$2,$4,$5,"Pacbio,Illumina","Pacbio"}}' > {output.sr_filter_anno}
        bcftools view --threads {threads} -H {input.lr_select_vcf} | awk -v OFS='\\t' '{{print$1,$2,$4,$5,"Pacbio,Illumina","Illumina"}}' > {output.lr_select_anno}
        bcftools view --threads {threads} -H {input.sr_spec_vcf} | awk -v OFS='\\t' '{{print$1,$2,$4,$5,"Illumina","Illumina"}}' > {output.sr_spec_anno}
        bcftools view --threads {threads} -H {input.lr_spec_vcf} | awk -v OFS='\\t' '{{print$1,$2,$4,$5,"Illumina","Illumina"}}' > {output.lr_spec_anno}
        bcftools view --threads {threads} -H {input.sr_callset_indel_vcf} | awk -v OFS='\\t' '{{print$1,$2,$4,$5,"Pacbio","Pacbio"}}' > {output.sr_indel_anno}

        cat {output.sr_filter_anno} {output.lr_select_anno} {output.sr_spec_anno} {output.lr_spec_anno} {output.sr_indel_anno} | \
            sort -k 1,1V -k 2,2n | bgzip -@ {threads} -c > {output.consensus_anno}
        tabix -s1 -b2 -e2 {output.consensus_anno}
        """



rule final_merge:
    input:
        sr_filter_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_shared.sub_filter.vcf.gz",
        lr_select_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_shared.sub_select.vcf.gz",
        sr_spec_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.srs_specific.vcf",
        lr_spec_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.lrs_specific.vcf",
        sr_callset_indel_vcf = "c3_merge_snv/callset/CKCG.analysis_set.call_set.srs.indel.vcf.gz",
        consensus_anno = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.anno.gz"
    output:
        consensus_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.vcf.gz"
    threads: 16

    shell:
        """
        echo -e '##INFO=<ID=CallsetPlatform,Number=.,Type=String,Description="Platforms of callsets that called this genotype">\\n##INFO=<ID=GenotypePlatform,Number=.,Type=String,Description="Platforms of callsets that use as the final genotype">' > c3_merge_snv/callset/srs_lrs_compare/header2

        

        
        bcftools concat -a --threads {threads} \
            {input.sr_filter_vcf} {input.lr_select_vcf} {input.sr_spec_vcf} {input.lr_spec_vcf} {input.sr_callset_indel_vcf} | \
            bcftools annotate --threads {threads} -h c3_merge_snv/callset/srs_lrs_compare/header2 \
            -a {input.consensus_anno} \
            -c CHROM,POS,REF,ALT,INFO/CallsetPlatform,INFO/GenotypePlatform \
            -x ID,INFO,FORMAT/DS,FORMAT/GP \
            -Oz \
            -o {output.consensus_vcf}
        tabix -f {output.consensus_vcf}
        """

rule generate_final_sr_scaffold:
    input:
        consensus_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.vcf.gz",
        concat_shapeit4_filter_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.maf0.01.vcf.gz"
        
    output:
        sr_scaffold_info = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.srs_scaffold.info.gz",
        sr_scaffold_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.srs_scaffold.vcf.gz"
    shell:
        """
        bcftools view --threads 16 -v snps -i "GenotypePlatform=='Illumina'" {input.consensus_vcf} | bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\tSCAFFOLD\\n" | bgzip -c > {output.sr_scaffold_info}

        echo '##FILTER=<ID=SCAFFOLD,Description="The variant in the short-read phasing scaffold">' > c3_merge_snv/callset/srs_lrs_compare/header3
        
        bcftools annotate --threads 16 -h c3_merge_snv/callset/srs_lrs_compare/header3 -a {output.sr_scaffold_info} -c CHROM,POS,REF,ALT,FILTER {input.concat_shapeit4_filter_vcf} | \
            bcftools view --threads 16 -f SCAFFOLD -o {output.sr_scaffold_vcf}
        """

#bcftools view --threads 16 -v snps -i "GenotypePlatform=='Illumina'" CKCG.analysis_set.call_set.consensus.vcf.gz | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\tSCAFFOLD\n" > CKCG.analysis_set.call_set.consensus.srs_scaffold.info
#bgzip CKCG.analysis_set.call_set.consensus.srs_scaffold.info

#bcftools view --threads 16 -v snps -i "CallsetPlatform=='Pacbio' && CallsetPlatform=='Illumina'" CKCG.analysis_set.call_set.consensus.vcf.gz -o CKCG.analysis_set.call_set.consensus.shared.vcf.gz
#tabix CKCG.analysis_set.call_set.consensus.shared.vcf.gz


generate_sample_consensus_vcf:
    input:
        consensus_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.vcf.gz"
    output:
        sample_consensus_vcf = "c3_merge_snv/samples/{sample}/{sample}.consensus.vcf",
        sample_consensus_miss_vcf = "c3_merge_snv/samples/{sample}/{sample}.consensus.miss.vcf.gz"
    threads: 
    shell:
        """
        bcftools view -s {wildcards.sample} --threads {threads} {input.consensus_vcf} | \
            bcftools view --threads {threads} \
            -e 'GT=="0/0" || GT=="0|0" || GT=="0"' | \
            awk -v OFS='\\t' '{{if(substr($1,1,1)=="#") print$0; else{{if(length($10)==1) {{gt=$10"/"$10; print $1,$2,$3,$4,$5,$6,$7,$8,$9,gt }} else print$0}} }}' \
            > {output.sample_consensus_vcf}

        bcftools view --threads {threads} -i "GT=='./.' || GT=='.|.'" {output.sample_consensus_vcf} -o {output.sample_consensus_miss_vcf}
        tabix {output.sample_consensus_miss_vcf}
        """















