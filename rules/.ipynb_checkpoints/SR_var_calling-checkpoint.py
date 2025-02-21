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
        #20G
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
        gvcf_map = "c1_call_sr_snv/sr_gvcf/{config['prefix']}.analysis_set.gvcf.map"
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
#         interval_file = "/storage/yangjianLab/wangyifei/project/01.{config['prefix']}/10.SR_variant/CHM13/03.merge/analysis_set/CHM13.20mb.interval"
#     output:
#         touch("c1_call_sr_snv/interval_vcf/.intervals_generated")
#     run:
#         with open(input.interval_file) as f:
#             intervals_list = [line.split("\t")[0] for line in f]
#         config['intervals'] = intervals_list
        
        
rule GenomicsDB_GenotypeGVCFs_interval:
    input:
        gvcf_map = "c1_call_sr_snv/sr_gvcf/{config['prefix']}.analysis_set.gvcf.map" 
    output:
        vcf = "c1_call_sr_snv/interval_vcf/{interval}.raw_variant.vcf.gz"
    log:
        "logs/GenomicsDB_GenotypeGVCFs_interval/{interval}.log"
    resources:
        #20G
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
        merged_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.raw_variant.analysis_set.vcf.gz"
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
        vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.raw_variant.analysis_set.vcf.gz",
        ref = config['reference']['CHM13'],
        hapmap = config['GATK_Resource']['hapmap'],
        omni = config['GATK_Resource']['omni'],
        _1000G = config['GATK_Resource']['1000G'],
        dbsnp = config['GATK_Resource']['dbsnp'],
        known_indel = config['GATK_Resource']['known_indel']
    output:
        snp_recalibrated_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.snp_recalibrated.analysis_set.vcf.gz"
    log:
        "logs/merged_vcf_snp_VQSR/merged_vcf_snp_VQSR.log"
    resources:
        # 100G-200G
        max_mem_mb = 200000,
        min_mem_mb = 100000    
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
        vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.snp_recalibrated.analysis_set.vcf.gz",
        ref = config['reference']['CHM13'],
        mills = config['GATK_Resource']['mills'],
        axiomPoly = config['GATK_Resource']['axiomPoly'],
        dbsnp = config['GATK_Resource']['dbsnp']
    output:
        snp_indel_recalibrated_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.analysis_set.vcf.gz"
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

#snp_indel_recalibrated_filter_auto_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.auto_chr.vcf.gz"
#snp_indel_recalibrated_filter_chrX_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.chrX.vcf.gz"
rule gatk_vcf_filter:
    input:
        snp_indel_recalibrated_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.analysis_set.vcf.gz",
        ref = config['reference']['CHM13']
    output:
        snp_indel_recalibrated_filter_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.analysis_set.biallelic.vcf.gz",
        snp_indel_recalibrated_filter_hwe_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.vcf.gz"
    params:
        variant_samples_number = 1034
    shell:
        """
        bcftools view --threads 16 -f PASS {input.snp_indel_recalibrated_vcf} | bcftools norm --threads 16 -m -any -f {input.ref} | bcftools view --threads 16 -e "ALT=='*'" | bcftools plugin fill-tags --threads 16 | bcftools view --threads 16 -i 'HWE>=1e-6 || ExcHet>=0.1' -o {output.snp_indel_recalibrated_filter_vcf}
        
        bcftools view --threads 16 -i "HWE>=1e-6 && NS>1034" {output.snp_indel_recalibrated_filter_vcf} -o {output.snp_indel_recalibrated_filter_hwe_vcf}

        #bcftools view --threads 16 -e "CHR=='chrX' | CHR=='chrY' | CHR=='chrM'" {output.snp_indel_recalibrated_filter_vcf} -o {output.snp_indel_recalibrated_filter_auto_vcf}

        #bcftools view --threads 16 -i "CHR=='chrX'" {output.snp_indel_recalibrated_filter_vcf} -o {output.snp_indel_recalibrated_filter_chrX_vcf}
        """



rule male_whatshap_phase:
    input:
        snp_indel_recalibrated_filter_hwe_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.vcf.gz",
        zmw_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.zmw.pbmm2.bam",
        ref = config['reference']['CHM13']
    output:
        sample_auto_vcf = temp("c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.{sample}.auto_chr.vcf.gz"),
        sample_auto_whatshap_vcf = temp("c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.{sample}.auto_chr.vcf.gz"),
        sample_chrX_vcf = temp("c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.{sample}.chrX.vcf.gz"),
        sample_chrX_whatshap_vcf = temp("c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.{sample}.chrX.vcf.gz"),
        sample_whatshap_vcf = "c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.{sample}.vcf.gz"
    threads: 4
    resources:
        #20G
        mem_mb = 20000
    shell:
        """
        bcftools view --threads {threads} -e "CHR=='chrX' | CHR=='chrY' | CHR=='chrM'" -s {wildcards.sample}-WGS -a {input.snp_indel_recalibrated_filter_hwe_vcf} | bcftools view --threads {threads} -e 'GT=="0/0" | GT=="0|0"' | bcftools norm --threads {threads} -m +any -f {input.ref} | whatshap unphase - | bgzip -@ {threads} -c > {output.sample_auto_vcf}
            
        whatshap phase --ignore-read-groups --reference {input.ref} --output {output.sample_auto_whatshap_vcf} {output.sample_auto_vcf} {input.zmw_bam}

        tabix -f {output.sample_auto_whatshap_vcf}

        bcftools view --threads {threads} -s {wildcards.sample}-WGS -a {input.snp_indel_recalibrated_filter_hwe_vcf} chrX | bcftools view --threads {threads} -e 'GT=="0/0" || GT=="0|0" || GT=="0"' | awk -v OFS='\\t' '{{if(substr($1,1,1)=="#") print$0; else{{if($1=="chrX" && $2>=2394411 && $2<=153925834) {{for(i=1;i<=9;i++) printf $i"\\t"; split($10,info,":");info[1]=info[1]"/"info[1];for(j=1;j<=length(info)-1;j++) printf info[j]":";printf info[length(info)];print"" }} else print$0 }} }}' | bcftools norm --threads {threads} -m +any -f {input.ref} | whatshap unphase - | bgzip -@ {threads} -c > {output.sample_chrX_vcf}
        
        whatshap phase --ignore-read-groups --reference {input.ref} --output {output.sample_chrX_whatshap_vcf} {output.sample_chrX_vcf} {input.zmw_bam}

        tabix -f {output.sample_chrX_whatshap_vcf}

        bcftools concat --threads {threads} -a -O z -o {output.sample_whatshap_vcf} {output.sample_auto_whatshap_vcf} {output.sample_chrX_whatshap_vcf}

        tabix -f {output.sample_whatshap_vcf}
        
        fi
         """

rule female_whatshap_phase:
    input:
        snp_indel_recalibrated_filter_hwe_vcf = "c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.vcf.gz"
    output:
        sample_vcf = temp("c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.analysis_set.biallelic.{sample}.vcf.gz"),
        sample_whatshap_vcf = "c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.{sample}.vcf.gz"
    threads: 4
    resources:
        #20G
        mem_mb = 20000
    shell:
        """
        bcftools view --threads {threads} -e "CHR=='chrY' | CHR=='chrM'" -s {wildcards.sample}-WGS -a {input.snp_indel_recalibrated_filter_hwe_vcf} | bcftools view --threads {threads} -e 'GT=="0/0" | GT=="0|0"' | bcftools norm --threads {threads} -m +any -f {input.ref} | whatshap unphase - | bgzip -@ {threads} -c > {output.sample_vcf}

        whatshap phase --ignore-read-groups --reference {input.ref} --output {output.sample_whatshap_vcf} {output.sample_vcf} {input.zmw_bam}

        tabix -f {output.sample_whatshap_vcf}
        """

# based on the sex, select the rules.
def dynamic_whatshap_phase_input(wildcards):
    if get_sex(wildcards) == "male":
        return rules.male_whatshap_phase.output

    else:
        return rules.female_whatshap_phase.output

use rule $(dynamic_whatshap_phase_input) as dynamic_step with:
    output: "c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.{sample}.vcf.gz"

#TODO:maybe here add into config file?
rule merged_whatshap_vcf:
    input:
        expand("c1_call_sr_snv/whatshap/{sample}/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.{sample}.vcf.gz", sample=config['samples']),
        ref = config['reference']['CHM13']
    output:
        whatshap_vcf_list = "c1_call_sr_snv/whatshap/merge/vcf.list",
        merged_whatshap_vcf = "CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.vcf.gz",
        merged_whatshap_filter_vcf = "CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.analysis_set.biallelic.vcf.gz"
    shell:
        """
        ls c1_call_sr_snv/whatshap/*/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.analysis_set.biallelic.*.vcf.gz > {output.whatshap_vcf_list}
        bcftools merge --threads 16 -0 -m any \
            -l {output.whatshap_vcf_list} \
            -Oz \
            -o {output.merged_whatshap_vcf}

        bcftools view --threads 16 {output.merged_whatshap_vcf} | \
            bcftools norm --threads 16 -m -any -f {input.ref} | \
            bcftools plugin fill-tags --threads 16 | \
            bcftools view --threads 16 -v snps -e "N_PASS(GT=='0|1' || GT=='1|0')==0 && (AN-AC<=1 || AC<=1)" \
            -o {output.merged_whatshap_filter_vcf}
        """

rule whatshap_vcf_shapeit4_topmed_scaffold:
    input:
        genetic_map = "/storage/yangjianLab/wangyifei/resource/genetic_map/chm13_liftover/shapeit4/genetic_map_chm13_{chr}.reform.txt",
        topmed_east_asian = "/storage/yangjianLab/wangyifei/resource/TOPMed/east_asian/merge_vcf/{chr}/TOPMed_WGS_freeze.8.east_asian.merge.{chr}.filter.snps.liftover_chm13.shapeit4.vcf.gz",
        merged_whatshap_filter_vcf = "CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.analysis_set.biallelic.vcf.gz"
    output:
        topmed_scaffold_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.scaffold.shapeit4.analysis_set.biallelic.{chr}.vcf.gz"
    threads: 16
    resources:
        #150G
        mem_mb = 150000
    shell:
        """
        shapeit4 \
            --input {input.merged_whatshap_filter_vcf} \
            --map {input.genetic_map} \
            --region {wildcards.chr} --pbwt-depth 8 -T {threads} --sequencing --use-PS 0.0001 \
            --reference {input.topmed_east_asian} \
            --out {output.topmed_scaffold_vcf}
        
        tabix -f {output.topmed_scaffold_vcf}
        """
rule whatshap_vcf_shapeit4:
    input:
        genetic_map = "/storage/yangjianLab/wangyifei/resource/genetic_map/chm13_liftover/shapeit4/genetic_map_chm13_{chr}.reform.txt",
        merged_whatshap_filter_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.analysis_set.biallelic.vcf.gz",
        topmed_scaffold_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.scaffold.shapeit4.analysis_set.biallelic.{chr}.vcf.gz"
    output:
        shapeit4_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.{chr}.vcf.gz"
    threads: 16
    resources:
        #150G
        mem_mb = 150000
    shell:
        """
        shapeit4 \
            --input {input.merged_whatshap_filter_vcf} \
            --map {input.genetic_map} \
            --region {wildcards.chr} --pbwt-depth 8 -T {threads} --sequencing --use-PS 0.0001 \
            --scaffold {input.topmed_scaffold_vcf} \
            --out {output.shapeit4_vcf}
        
        tabix -f {output.shapeit4_vcf}
        """

concat -f vcf.list --threads 16 -O z -o CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.vcf.gz

#chr1-chr22,chrX
rule shapeit4_vcf_concat:
    input:
        shapeit4_vcfs = expand("c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.{chr}.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        shapeit4_vcf_list = "c1_call_sr_snv/shapeit4/vcf.list",
        concat_shapeit4_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.vcf.gz",
        concat_shapeit4_filter_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.maf0.01.vcf.gz"
    threads: 16
    resources:
        mem_mb = 30000
    shell:
        """
        ls c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.chr*.vcf.gz > {output.shapeit4_vcf_list}
        
        bcftools concat -f {output.shapeit4_vcf_list} \
            --threads 16 \
            -Oz \
            -o {output.concat_shapeit4_vcf}

        bcftools view --threads 16 \
            -i 'AF>=0.01 && AF<=0.99' \
            -o {output.concat_shapeit4_filter_vcf} \
            {output.concat_shapeit4_vcf}
        """

# change the name into "xxx-CLR" to be used in long read whatshap.
rule rename_samples_concat_shapeit4_vcf:
    input:
        concat_shapeit4_filter_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.maf0.01.vcf.gz"
    output:
        rename_concat_shapeit4_vcf = "c1_call_sr_snv/shapeit4/CKCG.gatk.variant_recalibrated.filter.hwe_missing_filter.whatshap.unphase_singleton_filter.topmed_eas.shapeit4.analysis_set.biallelic.maf0.01.rename.vcf.gz"
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        """
        zcat {input.concat_shapeit4_filter_vcf} | sed '/^#CHROM/s/WGS/CLR/g' | bgzip -c > {output.rename_concat_shapeit4_vcf}
        """
rule xxx:
    input:
    output:
    threads:
    resources:
        mem_mb = 
    shell:
        """
        """
rule xxx:
    input:
    output:
    threads:
    resources:
        mem_mb = 
    shell:
        """
        """
rule xxx:
    input:
    output:
    threads:
    resources:
        mem_mb = 
    shell:
        """
        """






