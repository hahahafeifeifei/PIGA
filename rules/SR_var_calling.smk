# rule call_sr_snv:
#     input:
#         expand("c1_call_sr_snv", sample=config['samples'])
        
def get_sr_bwa_map_input_fastqs(wildcards):
    return config["sr_fastqs"][wildcards.sample]

rule sr_bwa_map:
    input:
        fastqs = get_sr_bwa_map_input_fastqs,
        ref = config['reference']['CHM13'] 
    output:
        bam = "c1_call_sr_snv/sr_mapping/{sample}/{sample}.srt.bam",
        bai = "c1_call_sr_snv/sr_mapping/{sample}/{sample}.srt.bam.bai"
    log:
        "logs/sr_bwa_map/{sample}.log" 
    threads: 4
    shell:
        """
        (bwa mem -t {threads} -Y -R '@RG\tID:'{wildcards.sample}'\tSM:'{wildcards.sample}'\tPL:ILLUMINA' {input.ref} {input.fastqs[0]} {input.fastqs[1]} | \
        samtools view -uSh -@ {threads} - | \
        samtools sort -O bam -@ {threads} -o {output.bam} - ) 2> {log}
        samtools index -@ {threads} {output.bam}
        """

#TODO: also consider the unpaired R1 fastq file and R2 fastq file.  

rule sr_bam_dedup:
    input:
        bam = "c1_call_sr_snv/sr_mapping/{sample}/{sample}.srt.bam",
        ref = config['reference']['CHM13']
    output:
        md_bam = "c1_call_sr_snv/sr_dedup/{sample}/{sample}.srt.dedup.bam",
        md_bai = "c1_call_sr_snv/sr_dedup/{sample}/{sample}.srt.dedup.bam.bai",
        md_metric = "c1_call_sr_snv/sr_dedup/{sample}/{sample}.srt.dedup.metrics"
    log:
        "logs/sr_bam_dedup/{sample}.log" 
    threads: 4
    shell:
        """
        gatk MarkDuplicates -R {input.ref} -I {input.bam} -O {output.md_bam} -M {output.md_metric} 2> {log}
        samtools index -@ {threads} {output.md_bam}
        """

# def get_sr_bam_BQSR_GATK_Resource(wildcards, config, rule):
#     return config["GATK_Resource"][rule.input.ref]

rule sr_bam_BQSR:
    input:
        md_bam = "c1_call_sr_snv/sr_dedup/{sample}/{sample}.srt.dedup.bam",
        md_metric = "c1_call_sr_snv/sr_dedup/{sample}/{sample}.srt.dedup.metrics",
        Mills_1000G = config['GATK_Resource']['mills'],
        dbsnp = config['GATK_Resource']['dbsnp'],
        known_indel = config['GATK_Resource']['known_indel'],
        ref = config['reference']['CHM13']
    output:
        bqsr_bam = "c1_call_sr_snv/sr_bqsr/{sample}/{sample}.bqsr.bam"
    log:
        "logs/sr_bam_BQSR/{sample}.log"    
    shell:
        """
        gatk BaseRecalibrator \
        -I {input.md_bam} \
        -R {input.ref} \
        --known-sites {input.known_indel} \
        --known-sites {input.Mills_1000G} \
        --known-sites {input.dbsnp} \
        -O c1_call_sr_snv/sr_bqsr/{wildcards.sample}/{wildcards.sample}.recal.table \
        2> {log}
        
        gatk ApplyBQSR \
        -R {input.ref} \
        -I {input.md_bam} \
        -O {output.bqsr_bam} \
        c1_call_sr_snv/sr_bqsr/{wildcards.sample}/{wildcards.sample}.recal.table \
        2>> {log}
        """


#TODO: consider the variant calling of chrX, chrY and also chrM.
rule HaplotypeCaller_autosomes:
    input:
        bam = "c1_call_sr_snv/sr_bqsr/{sample}/{sample}.bqsr.bam",
        ref = config['reference']['CHM13'],
        dbsnp = config['GATK_Resource']['dbsnp']
    output:
        gvcf = "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.g.vcf.gz"
    log:
        "logs/HaplotypeCaller_autosomes/{sample}.log"
    params:
        ploidy = 2,
        intervals = " ".join([f"-L chr{i}" for i in range(1, 23)])
    threads: 2
    resources:
        mem_mb = 20000
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_mb}M" HaplotypeCaller \
            -ploidy {params.ploidy} {params.intervals} \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -D {input.dbsnp} \
            -ERC GVCF \
            2> {log}
        tabix -f {output.gvcf}
        """
# rule merge_vcfs:
#     input:
#         ref = config['reference']['CHM13']
#         autosomes = "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.auto_chr.g.vcf.gz",
#         chrX_PAR = "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.chrX_PAR.g.vcf.gz",
#         chrX_nonPAR = "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.chrX_nonPAR.g.vcf.gz",
#         chrY = "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.chrY.g.vcf.gz",
#         chrM = "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.chrM.g.vcf.gz"
#     output:
#         vcf = "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.g.vcf.gz"
#     log:
#         "logs/HaplotypeCaller_autosomes/{sample}.log"
#     threads: 2
#     resources:
#         mem_mb = 20000
#     shell:
#         """
#         gatk --java-options "-Xmx{resources.mem_mb}M" MergeVcfs \
#             -R {input.ref} \
#             -I {input.autosomes} \
#             -I {input.chrX_PAR} \
#             -I {input.chrX_nonPAR} \
#             -I {input.chrY} \
#             -I {input.chrM} \
#             -O {output.vcf} \
#             2> {log}
#         """        

rule construct_gvcf_map:
    input:
        gvcf = expand("c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.g.vcf.gz",sample = config['samples'])
    output:
        gvcf_map = "c1_call_sr_snv/sr_gvcf/CKCG.analysis_set.gvcf.map"
    threads: 1
    shell:
        """
        ls c1_call_sr_snv/sr_gvcf/*/*.gatk.g.vcf.gz > {output.gvcf_map}
        """

# rule generate_intervals:
#     input:
#         []
#     output:
#         "CHM13.20mb.interval"
#     shell:
#         bedtools makewindows -g {input} -w 20000000 > {output}

# rule add_intervals_to_config:
#     input:
#         interval_file = "/storage/yangjianLab/wangyifei/project/01.CKCG/10.SR_variant/CHM13/03.merge/analysis_set/CHM13.20mb.interval"
#     output:
#         touch("c1_call_sr_snv/interval_vcf/.intervals_generated")
#     run:
#         with open(input.interval_file) as f:
#             intervals_list = [line.split("\t")[0] for line in f]
#         config['intervals'] = intervals_list
        
        
rule GenomicsDB_GenotypeGVCFs_interval:
    input:
        gvcf_map = "c1_call_sr_snv/sr_gvcf/CKCG.analysis_set.gvcf.map" 
    output:
        vcf = "c1_call_sr_snv/interval_vcf/{interval}.raw_variant.vcf.gz"
    log:
        "logs/GenomicsDB_GenotypeGVCFs_interval/{interval}.log"
    resources:
        mem_mb = 20000
    threads: 2
    params:
        ref_name = os.path.splitext(config["reference"]['CHM13'])[0],
        batch_size = 50
    shell:
        """
        rm -rf /data/{wildcards.interval}_tmp
        mkdir /data/{wildcards.interval}_tmp
        cp {params.ref_name}.fasta /data/{wildcards.interval}_tmp
        cp {params.ref_name}.fasta.fai /data/{wildcards.interval}_tmp
        cp {params.ref_name}.dict /data/{wildcards.interval}_tmp
        gatk --java-options "-Xmx{resources.mem_mb}m" GenomicsDBImport \
            --sample-name-map {input.gvcf_map} \
            --genomicsdb-workspace-path /data/{wildcards.interval}_db \
            -R /data/{wildcards.interval}_tmp/chm13v2.0_maskedY_rCRS.fasta \
            --batch-size {params.batch_size} \
            --reader-threads {threads} \
            -L {wildcards.interval} \
            --overwrite-existing-genomicsdb-workspace true \
            --tmp-dir /data/{wildcards.interval}_tmp 2> {log}
        rm -rf /data/{wildcards.interval}_tmp
        
        gatk --java-options "-Xmx20g" GenotypeGVCFs \
            -AX ExcessHet -AX InbreedingCoeff \
            -V gendb:///data/{wildcards.interval}_db \
            -R {params.ref_name}.fasta \
            -O {output.vcf} 2>> {log}
        rm -rf /data/{wildcards.interval}_db
        """


rule merge_intervals:
    input:
        vcfs = expand("c1_call_sr_snv/interval_vcf/{interval}.raw_variant.vcf.gz",interval = config['intervals'])
    output:
        merged_vcf = "c1_call_sr_snv/merged_vcf/CKCG.gatk.raw_variant.analysis_set.vcf.gz"
    log:
        "logs/merge_intervals/merge_intervals.log"
    shell:
        """
        ls c1_call_sr_snv/interval_vcf/*.raw_variant.vcf.gz > c1_call_sr_snv/merged_vcf/interval.vcf.list
        gatk MergeVcfs \
            -I c1_call_sr_snv/merged_vcf/interval.vcf.list \
            -O {output.merged_vcf} 2> {log}
        """

rule merged_vcf_snp_VQSR:
    input:
        vcf = "c1_call_sr_snv/merged_vcf/CKCG.gatk.raw_variant.analysis_set.vcf.gz",
        ref = config['reference']['CHM13'],
        hapmap = config['GATK_Resource']['hapmap'],
        omni = config['GATK_Resource']['omni'],
        _1000G = config['GATK_Resource']['1000G'],
        dbsnp = config['GATK_Resource']['dbsnp'],
        known_indel = config['GATK_Resource']['known_indel']
    output:
        snp_recalibrated_vcf = "c1_call_sr_snv/merged_vcf/CKCG.gatk.snp_recalibrated.analysis_set.vcf.gz"
    log:
        "logs/merged_vcf_snp_VQSR/merged_vcf_snp_VQSR.log"
    resources:
        max_mem_gb = 200,
        min_mem_gb = 100    
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}g -Xms{resources.min_mem_gb}g" \
            VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
            --resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input._1000G} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 {input.dbsnp} \
            -O c1_call_sr_snv/merged_vcf/snps.recal \
            --tranches-file c1_call_sr_snv/merged_vcf/snps.tranches \
            2> {log}
        
        gatk --java-options "-Xmx{resources.max_mem_gb}g -Xms{resources.min_mem_gb}g" \
            ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --recal-file c1_call_sr_snv/merged_vcf/snps.recal \
            --tranches-file c1_call_sr_snv/merged_vcf/snps.tranches \
            --truth-sensitivity-filter-level 99.5 \
            --create-output-variant-index true \
            -mode SNP \
            -O {output.snp_recalibrated_vcf} \
            2>> {log}
        """

rule merged_vcf_indel_VQSR:
    input:
        vcf = "c1_call_sr_snv/merged_vcf/CKCG.gatk.snp_recalibrated.analysis_set.vcf.gz",
        ref = config['reference']['CHM13'],
        mills = config['GATK_Resource']['mills'],
        axiomPoly = config['GATK_Resource']['axiomPoly'],
        dbsnp = config['GATK_Resource']['dbsnp']
    output:
        snp_indel_recalibrated_vcf = "c1_call_sr_snv/merged_vcf/CKCG.gatk.variant_recalibrated.analysis_set.vcf.gz"
    log:
        "logs/merged_vcf_indel_VQSR/merged_vcf_indel_VQSR.log"
    resources:
        max_mem_gb = 200,
        min_mem_gb = 100    
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}g -Xms{resources.min_mem_gb}g" \
            VariantRecalibrator \
            -V {input.vcf} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            -mode INDEL \
            --max-gaussians 4 \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \
            --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 {input.axiomPoly} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
            -O c1_call_sr_snv/merged_vcf/indels.recal \
            --tranches-file c1_call_sr_snv/merged_vcf/indels.tranches \
            2> {log}
        
        gatk --java-options "-Xmx{resources.max_mem_gb}g -Xms{resources.min_mem_gb}g" \
            ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --recal-file c1_call_sr_snv/merged_vcf/indels.recal \
            --tranches-file c1_call_sr_snv/merged_vcf/indels.tranches \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode INDEL \
            -O {output.snp_indel_recalibrated_vcf} \
            2>> {log}
        """
