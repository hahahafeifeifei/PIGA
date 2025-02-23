def get_sr_fastqs(wildcards)):
    return config["sr_fastqs"][wildcards.sample]

def get_sr_unmapped_fastqs(wildcards):
    return config["sr_unmapped_fastqs"][wildcards.sample]

def get_zmw_pbmm2_map_input_fastqs(wildcards):
    return config["lr_zmw_fastqs"][wildcards.sample]        

def get_hifi_pbmm2_map_input_fastqs(wildcards):
    return config["lr_hifi_fastqs"][wildcards.sample]
    
def concat_final_phase_vcf_sex_specific_chrlist(wildcards):
    sex = config['sex'][wildcards.sample]
    chr_list = [f'chr{i}' for i in range(1, 23)]
    if gender == 'female':
        return chr_list + ['chrX']
    return chr_list

def get_sex(wildcards):
    sex = config['sex'][wildcards.sample]
    return sex 

rule personal_ref_bwa_map:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        ngs_R1_fastq = get_sr_fastqs[0],
        ngs_R2_fastq = get_sr_fastqs[1],
        ngs_R1_unmapped_fastq = get_sr_unmapped_fastqs[0],
        ngs_R2_unmapped_fastq = get_sr_unmapped_fastqs[1]
    output:
        ngs_ref_01_bam = temp("c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.01.srt.bam"),
        ngs_ref_02_bam = temp("c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.02.srt.bam"),
        ngs_ref_03_bam = temp("c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.03.srt.bam"),
        ngs_ref_bam = temp("c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.srt.bam")
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        bwa mem -t {threads} -Y -R '@RG\\tID:'{wildcards.sample}'\\tSM:'{wildcards.sample}'\\tPL:ILLUMINA' {input.consensus_fasta} {input.ngs_R1_fastq} {input.ngs_R2_fastq} | samtools view -uSh -@ {threads} - | samtools sort -O bam -@ {threads} -o {output.ngs_ref_01_bam} -

        bwa mem -t {threads} -Y -R '@RG\\tID:'{wildcards.sample}'\\tSM:'{wildcards.sample}'\\tPL:ILLUMINA' {input.consensus_fasta} {input.ngs_R1_unmapped_fastq} | samtools view -uSh -@ {threads} - | samtools sort -O bam -@ {threads} -o {output.ngs_ref_02_bam} -

        bwa mem -t {threads} -Y -R '@RG\\tID:'{wildcards.sample}'\\tSM:'{wildcards.sample}'\\tPL:ILLUMINA' {input.consensus_fasta} {input.ngs_R2_unmapped_fastq} | samtools view -uSh -@ {threads} - | samtools sort -O bam -@ {threads} -o {output.ngs_ref_03_bam} -

        samtools merge -@ {threads} -O BAM -c -p -f {output.ngs_ref_bam} {output.ngs_ref_01_bam} {output.ngs_ref_02_bam} {output.ngs_ref_03_bam}
        samtools index {output.ngs_ref_bam}
        """

rule personal_ref_bam_dedup:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        ngs_ref_bam = "c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.srt.bam"
    output:
        ngs_dedup_ref_bam = temp("c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.srt.dedup.bam"),
        ngs_ref_bam_dup_matric = "c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.dedup.metrics"
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        ~/software/gatk-4.2.6.1/gatk MarkDuplicates \
            -R {input.consensus_fasta} \
            -I {input.ngs_ref_bam} \
            -O {output.ngs_dedup_ref_bam} \
            -M {output.ngs_ref_bam_dup_matric}
        samtools index {output.ngs_dedup_ref_bam}
        """

rule personal_ref_dv:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        ngs_dedup_ref_bam = "c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.srt.dedup.bam"
    output:
        personal_ref_dv_vcf = "c6_draft_assembly/personal_ref_dv/{sample}/{sample}.af_pangenome_polish.deepvariant.vcf.gz",
        personal_ref_dv_filter_vcf = "c6_draft_assembly/personal_ref_dv/{sample}/{sample}.af_pangenome_polish.deepvariant.filter.vcf.gz"
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        mkdir -p c6_draft_assembly/personal_ref_dv/{wildcards.sample}/tmp
        
        singularity exec -B -B $(pwd):/project \
            -B $(pwd)/c6_draft_assembly/personal_ref_dv/{wildcards.sample}/tmp:/tmp \
            -B ${{REF_DIR}}:${{REF_DIR}} \
            ~/software/deepvariant/deepvariant_v1.3_sandbox \
            /opt/deepvariant/bin/run_deepvariant \
            --num_shards {threads} \
            --model_type=PACBIO \
            --ref=/project/{input.consensus_fasta} \
            --reads=/project/{input.ngs_dedup_ref_bam} \
            --output_vcf=/project/{output.personal_ref_dv_vcf} \
            --sample_name {wildcards.sample}
        
        bcftools view -@ {threads} -f PASS -i "GQ>20" -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX {output.personal_ref_dv_vcf} | \
            bcftools norm -m -any -f {input.consensus_fasta} -o {output.personal_ref_dv_filter_vcf}
        tabix -f {output.personal_ref_dv_filter_vcf}
        """



rule prepare_sample_consensus_vcf:
    input:
        concat_consensus_whatshap_shapeit4_vcf = "c4_phase_snv/shapeit4/CKCG.consensus.whatshap.shapeit4.analysis_set.vcf.gz"
    output:
        sample_vcf = temp("c5_personal_ref/snv_liftover/{sample}/{sample}.shapeit4.vcf.gz")
    shell:
        """
        bcftools view --threads {threads} \
            -s {wildcards.sample} \
            {input.concat_consensus_whatshap_shapeit4_vcf} | \
            bcftools view --threads {threads} \
            -i "GT!='RR'" \
            -o {output.sample_vcf}

        tabix -f {output.sample_vcf}
        """






rule sample_origin_snv_liftover:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        chain = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.chain",
        sample_vcf = "c5_personal_ref/snv_liftover/{sample}/{sample}.shapeit4.vcf.gz"
    output:
        sample_assembly_vcf = "c6_draft_assembly/snv_liftover/{sample}/{sample}.assembly.shapeit4.vcf.gz"
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        gatk CreateSequenceDictionary -R {input.consensus_fasta}
        gatk LiftoverVcf \
            -I {input.sample_vcf} \
            -O {output.sample_assembly_vcf} \
            -C {input.chain} \
            --REJECT /dev/null \
            -R {input.consensus_fasta} \
            --RECOVER_SWAPPED_REF_ALT

        tabix -f {output.sample_assembly_vcf}
        """

#merge varaints which from deepvariant and varaints which from liftover.
rule merge_variants:
    input:
        personal_ref_dv_filter_vcf = "c6_draft_assembly/personal_ref_dv/{sample}/{sample}.af_pangenome_polish.deepvariant.filter.vcf.gz",
        sample_assembly_vcf = "c6_draft_assembly/snv_liftover/{sample}/{sample}.assembly.shapeit4.vcf.gz"
    output:
        sample_assembly_merge_vcf = temp("c6_draft_assembly/snv_merge/{sample}/{sample}.assembly.merge.vcf.gz")
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        bcftools merge \
            -m none \
            --force-sample \
            {input.personal_ref_dv_filter_vcf} \
            {input.sample_assembly_vcf} |\
            awk -v OFS='\\t' '{{if(substr($0,1,2)=="##") print$0; else{{if(substr($0,1,1)=="#") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10; else {{split($10,a,":");split($11,b,":"); if(b[1]=="./.")print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}} }} }}' |\
            bgzip -@ {threads} -c > {output.sample_assembly_merge_vcf}
        
        tabix -f {output.sample_assembly_merge_vcf}
        """

rule merfin_filter_merged_variants:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        sample_WGS_meryl = "c3_merge_snv/meryl/{sample}/{sample}-WGS.meryl/merylIndex",
        sample_assembly_merge_vcf = "c6_draft_assembly/snv_merge/{sample}/{sample}.assembly.merge.vcf.gz"
    output:
        sample_assembly_merge_filter_vcf = "c6_draft_assembly/snv_merge/{sample}/{sample}.assembly.merge.filter.vcf"
    threads: 16
    resources:
        #90G
        mem_mb = 90000
    params:
        mems = 90
    shell:
        """
        merfin -filter \
            -sequence {input.consensus_fasta} \
            -readmers c3_merge_snv/meryl/{wildcards.sample}/{wildcards.sample}-WGS.meryl \
            -vcf {input.sample_assembly_merge_vcf} \
            -output c6_draft_assembly/snv_merge/{wildcards.sample}/{wildcards.sample}.assembly.merge \
            -threads {threads} \
            -memory {params.mems}

        rm -rf c6_draft_assembly/snv_merge/{wildcards.sample}/CHM13.af_pangenome.{wildcards.sample}_polish.fasta.meryl
        """

rule personal_ref_zmw_pbmm2:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        zmw_fastq = get_zmw_pbmm2_map_input_fastqs
    output:
        zmw_ref_bam = temp("c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.assembly.zmw.pbmm2.bam")
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        pbmm2 align -j {threads} -J 4 \
            {input.consensus_fasta} \
            {input.zmw_fastq} \
            {output.zmw_ref_bam} \
            --sort
        """



#TODO: why phase a vcf by vcf itself?
rule sample_assembly_vcf_whatshap:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        sample_assembly_merge_filter_vcf = "c6_draft_assembly/snv_merge/{sample}/{sample}.assembly.merge.filter.vcf",
        zmw_ref_bam = "c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.assembly.zmw.pbmm2.bam"
    output:
        sample_assembly_merge_filter_whatshap_vcf = "c6_draft_assembly/whatshap/{sample}/{sample}.assembly.merge.whatshap_phase.vcf.gz",
        sample_phase_block_list = "c6_draft_assembly/whatshap/{sample}/{sample}.phase_block.list"
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        whatshap phase \
            --ignore-read-groups \
            --reference {input.consensus_fasta} \
            --output {output.sample_assembly_merge_filter_whatshap_vcf} \
            {input.sample_assembly_merge_filter_vcf} \
            {input.zmw_ref_bam} {input.sample_assembly_merge_filter_vcf}

        tabix -f {output.sample_assembly_merge_filter_whatshap_vcf}

        whatshap stats \
            --block-list {output.sample_phase_block_list} \
            {output.sample_assembly_merge_filter_whatshap_vcf}
        """

rule largest_haplotype_extract:
    input:
        sample_assembly_merge_filter_whatshap_vcf = "c6_draft_assembly/whatshap/{sample}/{sample}.assembly.merge.whatshap_phase.vcf.gz",
        sample_phase_block_list = "c6_draft_assembly/whatshap/{sample}/{sample}.phase_block.list"
    output:
        sample_chr_assembly_merge_final_phase_vcf = "c6_draft_assembly/whatshap/{sample}/{sample}.assembly.merge.final_phase.{chr}.vcf.gz"
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        phase_ps=$(awk -v OFS='\\t' -v chr={wildcards.chr} 'BEGIN{{ps=0;number=0}} {{if($2==chr && $6>=number) {{ps=$3;number=$6}} }} END{{print ps}}' {output.sample_phase_block_list})
        bcftools view --threads {threads} -i "PS==${{phase_ps}}" {input.sample_assembly_merge_filter_whatshap_vcf} {wildcards.chr} -o {output.sample_chr_assembly_merge_final_phase_vcf}
        tabix -f {output.sample_chr_assembly_merge_final_phase_vcf}
        """

rule concat_final_phase_vcf:
    input:
        sample_chr_assembly_merge_final_phase_vcfs = expand("c6_draft_assembly/whatshap/{sample}/{sample}.assembly.merge.final_phase.{chr}.vcf.gz", chr=concat_final_phase_vcf_sex_specific_chrlist, allow_missing=True)
    output:
        sample_assembly_merge_final_phase_vcf = "c6_draft_assembly/whatshap/{sample}/{sample}.assembly.merge.final_phase.vcf.gz"
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        bcftools concat --threads {threads} \
            c6_draft_assembly/whatshap/{wildcards.sample}/{wildcards.sample}.assembly.merge.final_phase.*.vcf.gz \
            -Oz -o {output.sample_assembly_merge_final_phase_vcf}

        tabix -f {output.sample_assembly_merge_final_phase_vcf}
        """

rule correct_switch_error:
    input:
        sample_assembly_merge_final_phase_vcf = "c6_draft_assembly/whatshap/{sample}/{sample}.assembly.merge.final_phase.vcf.gz",
        sample_assembly_vcf = "c6_draft_assembly/snv_liftover/{sample}/{sample}.assembly.shapeit4.vcf.gz"
    output:
        sample_switch_correct_vcf = "c6_draft_assembly/switch_err_correct/{sample}/{sample}.assembly.switch_correct.vcf.gz",
        sample_switch_correct_bed = "c6_draft_assembly/switch_err_correct/{sample}/{sample}.assembly.switch_correct.bed",
        sample_switch_correct_dv_vcf = temp("c6_draft_assembly/switch_err_correct/{sample}/{sample}.assembly.switch_correct.deepvariant.vcf.gz"),
        sample_switch_correct_shapeit4_vcf = temp("c6_draft_assembly/switch_err_correct/{sample}/{sample}.assembly.switch_correct.shapeit4.vcf.gz"),
        sample_switch_correct_final_vcf = temp("c6_draft_assembly/switch_err_correct/{sample}/{sample}.assembly.merge.final_phase.switch_correct.vcf.gz")
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        bcftools merge -m none --force-sample \
            {input.sample_assembly_merge_final_phase_vcf} \
            {input.sample_assembly_vcf} | \
            bcftools view -i "GT[0]='RA'" | \
            python3 ~/software/script/switch_correct.py - c6_draft_assembly/switch_err_correct/{wildcards.sample}/{wildcards.sample}.assembly.switch_correct.vcf {output.sample_switch_correct_bed}

        bgzip -c c6_draft_assembly/switch_err_correct/{wildcards.sample}/{wildcards.sample}.assembly.switch_correct.vcf
        tabix -f {output.sample_switch_correct_vcf}

        bcftools view -e "AF>=0" {output.sample_switch_correct_vcf} -T ^{output.sample_switch_correct_bed} -o {output.sample_switch_correct_dv_vcf}
        tabix -f {output.sample_switch_correct_dv_vcf}
        bcftools view -i "AF>=0" {output.sample_switch_correct_vcf} -o {output.sample_switch_correct_shapeit4_vcf}
        tabix -f {output.sample_switch_correct_shapeit4_vcf}
        
        bcftools concat -a {output.sample_switch_correct_shapeit4_vcf} {output.sample_switch_correct_dv_vcf} -Oz -o {output.sample_switch_correct_final_vcf}
        tabix -f {output.sample_switch_correct_final_vcf}
        """

###de novo assembly

#TODO: where to get pac here?
rule personal_ref_CLR_all_pbmm2:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        pac = "pac=/storage/yangjianLab/wangyifei/project/01.CKCG/03.CCS+CLR/sample.list"
    output:
        all_fofn = "c6_draft_assembly/result/{sample}/{sample}.CLR_all.fofn",
        all_ref_unsort_bam = "c6_draft_assembly/result/{sample}/{sample}.CLR_all.origin.bam",
        all_ref_bam = temp("c6_draft_assembly/result/{sample}/{sample}.CLR_all.bam")
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        grep {wildcards.sample} {input.pac} | cut -f 2 | tr " " "\n" > {output.all_fofn}
        pbmm2 align -j {threads} {input.consensus_fasta} {output.all_fofn} {output.all_ref_unsort_bam}

        samtools sort -@ {threads} -m 3G -T {wildcards.sample} -O bam -o {output.all_ref_bam} {output.all_ref_unsort_bam}
        samtools index -@ {threads} {output.all_ref_bam}
        """

rule personal_ref_hifi_pbmm2:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        hifi_fastq = get_hifi_pbmm2_map_input_fastqs
    output:
        hifi_ref_bam = temp("c6_draft_assembly/result/{sample}/{sample}.assembly.hifi.bam")
    threads:
    resources:
        #
        mem_mb = 
    shell:
        """
        pbmm2 align -j {threads} \
            -J 4 \
            {input.consensus_fasta} \
            {input.hifi_fastq} \
            {output.hifi_ref_bam} \
            --sort --preset CCS
        
        """

#PAR_region
rule liftover_chrX_par_region:
    input:
        chrX_par = "/storage/yangjianLab/wangyifei/resource/Reference/CHM13/chm13v2.0_chrX_PAR.bed",
        chain = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.chain"
    output:
        chrX_assembly_par = "c6_draft_assembly/result/{sample}/{sample}.assembly.chrX_PAR.bed"
    shell:
        """
        liftOver {input.chrX_par} {input.chain} {output.chrX_assembly_par} /dev/null
        """



rule phase_assembly:
    input:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta",
        ngs_dedup_ref_bam = "c6_draft_assembly/personal_ref_mapping/{sample}/{sample}.srt.dedup.bam",
        all_ref_bam = "c6_draft_assembly/result/{sample}/{sample}.CLR_all.bam",
        hifi_ref_bam = "c6_draft_assembly/result/{sample}/{sample}.assembly.hifi.bam",
        sample_switch_correct_final_vcf = "c6_draft_assembly/switch_err_correct/{sample}/{sample}.assembly.merge.final_phase.switch_correct.vcf.gz",
        chrX_assembly_par = "c6_draft_assembly/result/{sample}/{sample}.assembly.chrX_PAR.bed"
    output:
    params:
        sex = get_sex,
        tmp_dir = "c6_draft_assembly/result/{sample}/{sample}_tmp"
    threads: 28
    resources:
        #200G
        mem_mb = 200000 
    shell:
        """
        python3 ~/software/script/phase-assembly/phase_assembly.version2.py \
            --ngs-bam {input.ngs_dedup_ref_bam} \
            --subread-bam {input.all_ref_bam} \
            --hifi-bam {input.hifi_ref_bam} \
            -r {input.consensus_fasta} \
            -v {input.sample_switch_correct_final_vcf} \
            -p {wildcards.sample} \
            --par-bed {input.chrX_assembly_par} \
            -g {params.sex} \
            -o c6_draft_assembly/result/{wildcards.sample}/assembly \
            -T {params.tmp_dir} \
            -t {threads}
        """

rule draft_assembly_adaptor_mask:
    input:
        hap1_fa = "c6_draft_assembly/result/{sample}/assembly/{sample}.hap1.fasta",
        hap2_fa = "c6_draft_assembly/result/{sample}/assembly/{sample}.hap2.fasta",
        adaptor_fa = "/storage/yangjianLab/duanzhongqu/software/HiFiAdapterFilt-2.0.0/CLR/adaptor.fa"
    output:
        hap1_adaptor_bed = "c6_draft_assembly/result/{sample}/assembly/{sample}.hap1.adaptor.bed",
        hap2_adaptor_bed = "c6_draft_assembly/result/{sample}/assembly/{sample}.hap2.adaptor.bed",
        hap1_adaptor_masked_fa = "c6_draft_assembly/result/{sample}/assembly/{sample}.hap1.adaptor_masked.fasta",
        hap2_adaptor_masked_fa = "c6_draft_assembly/result/{sample}/assembly/{sample}.hap2.adaptor_masked.fasta"
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    shell:
        """
        minimap2 -t {threads} -cxsr -f5000 -N2000 -secondary=yes \
            --cs {input.hap1_fa} \
            {input.adaptor_fa} | \
            awk -v OFS='\\t' '{{print$6,$8,$9}}' | \
            bedtools sort -i - | \
            bedtools merge -i - \
            > {output.hap1_adaptor_bed}
        
        bedtools maskfasta -fi {input.hap1_fa} -bed {output.hap1_adaptor_bed} -fo {output.hap1_adaptor_masked_fa}

        minimap2 -t {threads} -cxsr -f5000 -N2000 -secondary=yes \
            --cs {input.hap2_fa} \
            {input.adaptor_fa} | \
            awk -v OFS='\\t' '{{print$6,$8,$9}}' | \
            bedtools sort -i - | \
            bedtools merge -i - \
            > {output.hap2_adaptor_bed}
        
        bedtools maskfasta -fi {input.hap2_fa} -bed {output.hap2_adaptor_bed} -fo {output.hap2_adaptor_masked_fa}
        """





# rule xxx:
#     input:
#     output:
#     threads:
#     resources:
#         #
#         mem_mb = 
#     shell:
#         """
#         """
