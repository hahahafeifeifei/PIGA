rule all_call_sr_snv:
    input:
        f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.vcf.gz",
        
rule sr_bwa_map:
    input:
        sr_fq1 = config['sr_fastqs'][0],
        sr_fq2 = config['sr_fastqs'][1],
        ref = config['reference']['CHM13']
    output:
        bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.bam",
        bai = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.bam.bai"
    threads: 4
    shell:
        """
        bwa mem -t {threads} -Y -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' {input.ref} {input.sr_fq1} {input.sr_fq2} | \
        samtools view -uSh -@ {threads} - | \
        samtools sort -O bam -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """


rule sr_bam_dedup:
    input:
        bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.bam",
        ref = config['reference']['CHM13']
    output:
        md_bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.dedup.bam",
        md_bai = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.dedup.bam.bai",
        md_metric = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.dedup.metrics"
    threads: 4
    shell:
        """
        gatk MarkDuplicates -R {input.ref} -I {input.bam} -O {output.md_bam} -M {output.md_metric}
        samtools index -@ {threads} {output.md_bam}
        """

rule sr_bam_BQSR:
    input:
        md_bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.dedup.bam",
        md_metric = "c1_call_sr_snv/sample_bam/{sample}/{sample}.srt.dedup.metrics",
        Mills_1000G = config['GATK_Resource']['mills'],
        known_indel = config['GATK_Resource']['known_indel'],
        ref = config['reference']['CHM13']
    output:
        bqsr_bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.bqsr.bam",
        recal_table = "c1_call_sr_snv/sample_bam/{sample}/{sample}.recal.table"
    shell:
        """
        gatk BaseRecalibrator \
        -I {input.md_bam} \
        -R {input.ref} \
        --known-sites {input.known_indel} \
        --known-sites {input.Mills_1000G} \
        -O {output.recal_table}
        
        gatk ApplyBQSR \
        -R {input.ref} \
        -I {input.md_bam} \
        -O {output.bqsr_bam} \
        --bqsr-recal-file {output.recal_table}
        """

rule HaplotypeCaller_autosomes:
    input:
        bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.bqsr.bam",
        ref = config['reference']['CHM13']
    output:
        auto_gvcf = temp("c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.auto_chr.g.vcf.gz")
    params:
        ploidy = 2,
        intervals = " ".join([f"-L chr{i}" for i in range(1, 23)])
    threads: 2
    resources:
        max_mem_gb = 20
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G" HaplotypeCaller \
            -ploidy {params.ploidy} {params.intervals} \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.auto_gvcf} \
            -ERC GVCF

        tabix -f {output.auto_gvcf}
        """

rule HaplotypeCaller_male_chrX:
    input:
        bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.bqsr.bam",
        ref = config['reference']['CHM13']
    output:
        chrX_PAR_gvcf = temp("c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrX_PAR.g.vcf.gz"),
        chrX_nonPAR_gvcf = temp("c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrX_nonPAR.g.vcf.gz")
    params:
        PAR_ploidy = 2,
        nonPAR_ploidy = 1
    threads: 2
    resources:
        max_mem_gb = 20
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G" HaplotypeCaller \
            -ploidy {params.PAR_ploidy} \
            -L chrX:1-2394410 -L chrX:153925835-154259566 \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.chrX_PAR_gvcf} \
            -ERC GVCF
        tabix -f {output.chrX_PAR_gvcf}
        
        gatk --java-options "-Xmx{resources.max_mem_gb}G" HaplotypeCaller \
            -ploidy {params.nonPAR_ploidy} \
            -L chrX:2394411-153925834 \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.chrX_nonPAR_gvcf} \
            -ERC GVCF
        tabix -f {output.chrX_nonPAR_gvcf}
        """

rule HaplotypeCaller_female_chrX:
    input:
        bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.bqsr.bam",
        ref = config['reference']['CHM13']
    output:
        chrX_gvcf = temp("c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrX.g.vcf.gz")
    params:
        ploidy = 2
    threads: 2
    resources:
        max_mem_gb = 20
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G" HaplotypeCaller \
            -ploidy {params.ploidy} \
            -L chrX \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.chrX_gvcf} \
            -ERC GVCF
        tabix -f {output.chrX_gvcf}
        """


rule HaplotypeCaller_male_chrY:
    input:
        bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.bqsr.bam",
        ref = config['reference']['CHM13']
    output:
        chrY_gvcf = temp("c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrY.g.vcf.gz")
    params:
        ploidy = 1
    threads: 2
    resources:
        max_mem_gb = 20
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G" HaplotypeCaller \
            -ploidy {params.ploidy} \
            -L chrY \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.chrY_gvcf} \
            -ERC GVCF
        tabix -f {output.chrY_gvcf}
        """

rule HaplotypeCaller_chrM:
    input:
        bam = "c1_call_sr_snv/sample_bam/{sample}/{sample}.bqsr.bam",
        ref = config['reference']['CHM13']
    output:
        chrM_gvcf = temp("c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrM.g.vcf.gz")
    params:
        ploidy = 1
    threads: 2
    resources:
        max_mem_gb = 20
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G" HaplotypeCaller \
            -ploidy {params.ploidy} \
            -L chrM \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.chrM_gvcf} \
            -ERC GVCF
        tabix -f {output.chrM_gvcf}
        """


rule male_merge_vcfs:
    input:
        ref = config['reference']['CHM13'],
        autosomes = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.auto_chr.g.vcf.gz",
        chrX_PAR = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrX_PAR.g.vcf.gz",
        chrX_nonPAR = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrX_nonPAR.g.vcf.gz",
        chrY = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrY.g.vcf.gz",
        chrM = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrM.g.vcf.gz"
    output:
        vcf = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.male.gatk.g.vcf.gz"
    threads: 2
    resources:
        max_mem_gb = 20
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G" MergeVcfs \
            -R {input.ref} \
            -I {input.autosomes} \
            -I {input.chrX_PAR} \
            -I {input.chrX_nonPAR} \
            -I {input.chrY} \
            -I {input.chrM} \
            -O {output.vcf}
        """        

rule female_merge_vcfs:
    input:
        ref = config['reference']['CHM13'],
        autosomes = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.auto_chr.g.vcf.gz",
        chrX = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrX.g.vcf.gz",
        chrM = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.chrM.g.vcf.gz"
    output:
        vcf = "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.female.gatk.g.vcf.gz"
    threads: 2
    resources:
        max_mem_gb = 20
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G" MergeVcfs \
            -R {input.ref} \
            -I {input.autosomes} \
            -I {input.chrX} \
            -I {input.chrM} \
            -O {output.vcf}
        """        

#based on the sex, select the rules.
def dynamic_merge_vcfs_input(wildcards):
    if get_sex(wildcards) == "male":
        return rules.male_merge_vcfs.output

    else:
        return rules.female_merge_vcfs.output

# use rule $(dynamic_merge_vcfs_input) as dynamic_step with:
#     output: "c1_call_sr_snv/sr_gvcf/{sample}/{sample}.gatk.g.vcf.gz"

rule merge_vcfs:
  input:
    dynamic_merge_vcfs_input
  output:
    "c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.g.vcf.gz"
  shell:
    """
    mv {input} {output}
    tabix -f {output}
    """


rule construct_gvcf_map:
    input:
        gvcf = expand("c1_call_sr_snv/sample_gvcf/{sample}/{sample}.gatk.g.vcf.gz",sample = samples_list)
    output:
        gvcf_map = f"c1_call_sr_snv/sample_gvcf/{config['prefix']}.gvcf.map"
    threads: 1
    shell:
        """
        ls {input.gvcf} > {output.gvcf_map}
        awk -F "/" '{{print $(NF-1)"\\t"$0}}' {output.gvcf_map} > c1_call_sr_snv/sample_gvcf/tmp.gvcf.map && mv c1_call_sr_snv/sample_gvcf/tmp.gvcf.map {output.gvcf_map}
        """

        
rule GenomicsDB_GenotypeGVCFs_interval:
    input:
        gvcf_map = f"c1_call_sr_snv/sample_gvcf/{config['prefix']}.gvcf.map",
        ref = config['reference']['CHM13']
    output:
        vcf = "c1_call_sr_snv/interval_vcf/{interval}.raw_variant.vcf.gz"
    resources:
        max_mem_gb = 20
    threads: 2
    params:
        interval_tmp = "c1_call_sr_snv/interval_vcf/{interval}_tmp",
        interval_db = "c1_call_sr_snv/interval_vcf/{interval}_db",
        batch_size = 50
    shell:
        """
        rm -rf {params.interval_tmp}
        mkdir {params.interval_tmp}
        gatk --java-options "-Xmx{resources.max_mem_gb}G" GenomicsDBImport \
            --sample-name-map {input.gvcf_map} \
            --genomicsdb-workspace-path {params.interval_db} \
            -R {input.ref} \
            --batch-size {params.batch_size} \
            --reader-threads {threads} \
            -L {wildcards.interval} \
            --overwrite-existing-genomicsdb-workspace true \
            --tmp-dir {params.interval_tmp}
        rm -rf {params.interval_tmp}
        
        gatk --java-options "-Xmx20g" GenotypeGVCFs \
            -AX ExcessHet -AX InbreedingCoeff \
            -V gendb://{params.interval_db} \
            -R {input.ref} \
            -O {output.vcf}
        rm -rf {params.interval_db}
        """

checkpoint generate_intervals:
    input:
        ref_fai = config['reference']['CHM13'] + ".fai"
    output:
        genome_size = "c1_call_sr_snv/CHM13.size",
        genome_interval = "c1_call_sr_snv/CHM13.20mb.interval"
    shell:
        """
        awk -v FS='\t' '{{print $1"\\t"$2}}' {input.ref_fai} > {output.genome_size}
        bedtools makewindows -g {output.genome_size} -w 20000000 | awk '{print $1":"$2+1"-"$3}' > {output.genome_interval}
        """

def get_intervals_list(wildcards):
    genome_interval_file = checkpoints.generate_intervals.get().output[1]
    with open(genome_interval_file) as f:
        intervals_list = [line.split("\t")[0].strip() for line in f]
    return intervals_list


rule merge_intervals:
    input:
        vcfs = expand("c1_call_sr_snv/interval_vcf/{interval}.raw_variant.vcf.gz",interval = get_intervals_list)
    output:
        vcf_list = "c1_call_sr_snv/merged_vcf/interval.vcf.list",
        merged_vcf = f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.raw_variant.vcf.gz"
    shell:
        """
        ls {input.vcfs} > {output.vcf_list}
        gatk MergeVcfs \
            -I {output.vcf_list} \
            -O {output.merged_vcf}
        """

rule merged_vcf_snp_VQSR:
    input:
        vcf = f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.raw_variant.vcf.gz",
        ref = config['reference']['CHM13'],
        hapmap = config['GATK_Resource']['hapmap'],
        omni = config['GATK_Resource']['omni'],
        _1000G = config['GATK_Resource']['1000G'],
        known_indel = config['GATK_Resource']['known_indel']
    output:
        snp_recalibrated_vcf = f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.snp_recalibrated.vcf.gz",
        snp_recal = "c1_call_sr_snv/merged_vcf/snps.recal",
        snp_tranches = "c1_call_sr_snv/merged_vcf/snps.tranches"
    resources:
        max_mem_gb = 60,
        min_mem_gb = 40    
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G -Xms{resources.min_mem_gb}G" \
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
            -O {output.snp_recal} \
            --tranches-file {output.snp_tranches} \
        
        gatk --java-options "-Xmx{resources.max_mem_gb}G -Xms{resources.min_mem_gb}G" \
            ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --recal-file {output.snp_recal} \
            --tranches-file {output.snp_tranches} \
            --truth-sensitivity-filter-level 99.5 \
            --create-output-variant-index true \
            -mode SNP \
            -O {output.snp_recalibrated_vcf} \
        """

rule merged_vcf_indel_VQSR:
    input:
        vcf = f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.snp_recalibrated.vcf.gz",
        ref = config['reference']['CHM13'],
        mills = config['GATK_Resource']['mills'],
        axiomPoly = config['GATK_Resource']['axiomPoly']
    output:
        snp_indel_recalibrated_vcf = f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.vcf.gz",
        indel_recal = "c1_call_sr_snv/merged_vcf/indels.recal",
        indel_tranches = "c1_call_sr_snv/merged_vcf/indels.tranches"
    resources:
        max_mem_gb = 60,
        min_mem_gb = 40    
    shell:
        """
        gatk --java-options "-Xmx{resources.max_mem_gb}G -Xms{resources.min_mem_gb}G" \
            VariantRecalibrator \
            -V {input.vcf} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            -mode INDEL \
            --max-gaussians 4 \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \
            --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 {input.axiomPoly} \
            -O {output.indel_recal} \
            --tranches-file {output.indel_tranches} \
        
        gatk --java-options "-Xmx{resources.max_mem_gb}G -Xms{resources.min_mem_gb}G" \
            ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --recal-file {output.indel_recal} \
            --tranches-file {output.indel_tranches} \
            --truth-sensitivity-filter-level 99.0 \
            --create-output-variant-index true \
            -mode INDEL \
            -O {output.snp_indel_recalibrated_vcf} \
        """

#snp_indel_recalibrated_filter_auto_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.auto_chr.vcf.gz"
#snp_indel_recalibrated_filter_chrX_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.chrX.vcf.gz"
rule gatk_vcf_filter:
    input:
        snp_indel_recalibrated_vcf = f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.vcf.gz",
        ref = config['reference']['CHM13']
    output:
        snp_indel_recalibrated_filter_vcf = f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.vcf.gz",
    params:
        missing = len(samples_list)
    threads: 16
    shell:
        """
        bcftools view --threads {threads} -f PASS {input.snp_indel_recalibrated_vcf} | \
        bcftools annotate -x QUAL,INFO,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PGT,FORMAT/PID,FORMAT/PL,FORMAT/PS | \
        bcftools norm --threads {threads} -m -any -f {input.ref} | \
        bcftools view --threads {threads} -e "ALT=='*'" | \
        bcftools plugin fill-tags --threads {threads} | \
        bcftools view --threads {threads} -v snps -e "(HWE<1e-6 && ExcHet<0.1) || F_MISSING>0.1 || AC<1" -o {output.snp_indel_recalibrated_filter_vcf}
        tabix {output.snp_indel_recalibrated_filter_vcf}
        """



