rule all_draft_assembly:
    input:
        expand("c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap1.rename.fasta", sample = samples_list),
        expand("c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap2.rename.fasta", sample = samples_list)

# Create a symblink to the consensus fasta.
rule personal_ref_create_dict:
    input:
        consensus_fasta = get_consensus_fasta_input
    output:
        consensus_fasta_dict = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.dict"
    resources:
        mem_mb = 50000
    run:
        if "consensus_fasta" in config:
            consensus_fasta = f"c5_personal_ref/sample_reference/{wildcards.sample}/{wildcards.sample}.personal_ref.fasta"
            os.symlink(config['consensus_fasta'], consensus_fasta)
        else:
            # Here input.consensus_fasta is actually c5_personal_ref/sample_reference/{wildcards.sample}/{wildcards.sample}.personal_ref.fasta
            consensus_fasta = input.consensus_fasta

        shell(f"""
              gatk CreateSequenceDictionary -R {consensus_fasta}
        """)

rule personal_ref_bwa_map:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        consensus_fasta_dict = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.dict",
        sr_fq1 = config['sr_fastqs'][0],
        sr_fq2 = config['sr_fastqs'][1]
    output:
        ngs_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.srt.bam"
    threads: 16
    resources:
        mem_mb = 64*1024
    shell:
        """
        bwa index {input.consensus_fasta}
        samtools faidx {input.consensus_fasta}
        bwa mem -t {threads} -Y -R '@RG\\tID:'{wildcards.sample}'\\tSM:'{wildcards.sample}'\\tPL:ILLUMINA' {input.consensus_fasta} {input.sr_fq1} {input.sr_fq2} | \
        samtools view -uSh -@ {threads} - | \
        samtools sort -O bam -@ {threads} -o {output.ngs_bam}
        samtools index -@ {threads} {output.ngs_bam}
        """

rule personal_ref_bam_dedup:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        ngs_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.srt.bam"
    output:
        ngs_dedup_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.srt.dedup.bam",
        ngs_dedup_matric = "c6_draft_assembly/sample_assembly/{sample}/{sample}.dedup.metrics"
    threads: 4
    resources:
        mem_mb = 30*1024
    shell:
        """
       gatk MarkDuplicates \
            -R {input.consensus_fasta} \
            -I {input.ngs_bam} \
            -O {output.ngs_dedup_bam} \
            -M {output.ngs_dedup_matric}
        samtools index {output.ngs_dedup_bam}
        """

rule sample_origin_snv_liftover:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        consensus_fasta_dict = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.dict",
        chain = get_chain_input,
        sample_vcf = get_sample_vcf_input
    output:
        sample_personal_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.shapeit.vcf.gz"
    threads: 1
    resources:
        mem_mb = 30*1024
    shell:
        """
        gatk LiftoverVcf \
            -I {input.sample_vcf} \
            -O {output.sample_personal_vcf} \
            -C {input.chain} \
            --REJECT /dev/null \
            -R {input.consensus_fasta} \
            --RECOVER_SWAPPED_REF_ALT

        tabix -f {output.sample_personal_vcf}
        """

rule personal_ref_dv:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        ngs_dedup_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.srt.dedup.bam"
    output:
        personal_dv_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.deepvariant.vcf.gz",
        personal_dv_filter_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.deepvariant.filter.vcf.gz"
    singularity:
        "scripts/call_lr_snv/deepvariant.sif"
    threads: 16
    resources:
        mem_mb = 64*1024
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
            --num_shards {threads} \
            --model_type=PACBIO \
            --ref={input.consensus_fasta} \
            --reads={input.ngs_dedup_bam} \
            --output_vcf={output.personal_dv_vcf} \
            --sample_name {wildcards.sample}
        
        bcftools view --threads {threads} -f PASS -i "GQ>20" -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX {output.personal_dv_vcf} | \
            bcftools norm -m -any -f {input.consensus_fasta} -Oz -o {output.personal_dv_filter_vcf}
        tabix -f {output.personal_dv_filter_vcf}
        """

#merge varaints which from deepvariant and varaints which from liftover.
rule merge_variants:
    input:
        personal_dv_filter_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.deepvariant.filter.vcf.gz",
        sample_personal_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.shapeit.vcf.gz"
    output:
        sample_personal_merge_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.vcf.gz"
    threads: 4
    resources:
        mem_mb = 30*1024
    shell:
        """
        bcftools merge \
            -m none --force-sample \
            {input.personal_dv_filter_vcf} \
            {input.sample_personal_vcf} |\
            awk -v OFS='\\t' '{{if(substr($0,1,2)=="##") print$0; else{{if(substr($0,1,1)=="#") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10; else {{split($10,a,":");split($11,b,":"); if(b[1]=="./.")print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}} }} }}' |\
            bgzip -@ {threads} -c > {output.sample_personal_merge_vcf}
        
        tabix -f {output.sample_personal_merge_vcf}
        """

rule personal_ref_zmw_pbmm2:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        lr_zmw_fastqs = config['lr_zmw_fastqs']
    output:
        zmw_ref_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.zmw.srt.bam"
    threads: 8
    resources:
        mem_mb = 64*1024
    shell:
        """
        pbmm2 align -j {threads} -J 1 \
            {input.consensus_fasta} \
            {input.lr_zmw_fastqs} \
            {output.zmw_ref_bam} \
            --sort
        """

rule sample_assembly_vcf_whatshap:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        sample_assembly_merge_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.vcf.gz",
        zmw_ref_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.zmw.srt.bam"
    output:
        sample_assembly_merge_whatshap_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.phase.vcf.gz",
        sample_phase_block_list = "c6_draft_assembly/sample_assembly/{sample}/{sample}.phase_block.list"
    threads: 4
    resources:
        mem_mb = 30*1024
    shell:
        """
        whatshap phase \
            --ignore-read-groups \
            --reference {input.consensus_fasta} \
            --output {output.sample_assembly_merge_whatshap_vcf} \
            {input.sample_assembly_merge_vcf} \
            {input.zmw_ref_bam} {input.sample_assembly_merge_vcf}

        tabix -f {output.sample_assembly_merge_whatshap_vcf}

        whatshap stats \
            --block-list {output.sample_phase_block_list} \
            {output.sample_assembly_merge_whatshap_vcf}
        """

rule haplotype_largest_block:
    input:
        sample_assembly_merge_whatshap_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.whatshap_phase.vcf.gz",
        sample_phase_block_list = "c6_draft_assembly/sample_assembly/{sample}/{sample}.phase_block.list"
    output:
        sample_chr_assembly_merge_phase_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.whatshap_phase.{chr}.vcf.gz"
    threads: 4
    resources:
        mem_mb = 30*1024
    shell:
        """
        phase_ps=$(awk -v OFS='\\t' -v chr={wildcards.chr} 'BEGIN{{ps=0;number=0}} {{if($2==chr && $6>=number) {{ps=$3;number=$6}} }} END{{print ps}}' {input.sample_phase_block_list})
        bcftools view --threads {threads} -i "PS==${{phase_ps}}" {input.sample_assembly_merge_whatshap_vcf} {wildcards.chr} -o {output.sample_chr_assembly_merge_phase_vcf}
        tabix -f {output.sample_chr_assembly_merge_phase_vcf}
        """

rule concat_final_phase_vcf:
    input:
        sample_chr_assembly_merge_phase_vcfs = expand("c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.whatshap_phase.{chr}.vcf.gz", chr=concat_final_phase_vcf_sex_specific_chrlist, allow_missing=True)
    output:
        sample_assembly_merge_phase_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.phase.vcf.gz"
    threads: 4
    resources: 
        mem_mb = 30*1024
    shell:
        """
        bcftools concat --threads {threads} {input.sample_chr_assembly_merge_phase_vcfs} \
            -Oz -o {output.sample_assembly_merge_phase_vcf}
        tabix -f {output.sample_assembly_merge_phase_vcf}
        """

rule correct_switch_error:
    input:
        sample_assembly_merge_phase_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.phase.vcf.gz",
        sample_assembly_shapeit_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.shapeit.vcf.gz"
    output:
        sample_switch_correct_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.switch_correct.vcf.gz",
        sample_switch_correct_bed = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.switch_correct.bed",
        sample_switch_correct_dv_vcf = temp("c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.switch_correct.deepvariant.vcf.gz"),
        sample_switch_correct_shapeit_vcf = temp("c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.switch_correct.shapeit.vcf.gz"),
        sample_switch_correct_final_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.phase.switch_correct.vcf.gz"
    params:
        sample_switch_correct_unzip_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.switch_correct.vcf"
    threads: 4
    resources:
        mem_mb = 30*1024
    shell:
        """
        bcftools merge -m none --force-sample \
            {input.sample_assembly_merge_phase_vcf} \
            {input.sample_assembly_shapeit_vcf} | \
            bcftools view -i "GT[0]='RA'" | \
            python3 scripts/draft_assembly/switch_correct.py - {params.sample_switch_correct_unzip_vcf} {output.sample_switch_correct_bed}

        bgzip -f {params.sample_switch_correct_unzip_vcf}
        tabix -f {output.sample_switch_correct_vcf}

        bcftools view -e "AF>=0" {output.sample_switch_correct_vcf} -T ^{output.sample_switch_correct_bed} -o {output.sample_switch_correct_dv_vcf}
        tabix -f {output.sample_switch_correct_dv_vcf}
        bcftools view -i "AF>=0" {output.sample_switch_correct_vcf} -o {output.sample_switch_correct_shapeit_vcf}
        tabix -f {output.sample_switch_correct_shapeit_vcf}
        
        bcftools concat -a {output.sample_switch_correct_shapeit_vcf} {output.sample_switch_correct_dv_vcf} -Oz -o {output.sample_switch_correct_final_vcf}
        tabix -f {output.sample_switch_correct_final_vcf}
        """


rule personal_ref_subreads_pbmm2:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        lr_subreads_bam = config['lr_subreads_bam']
    output:
        subreads_ref_unsort_bam = temp("c6_draft_assembly/sample_assembly/{sample}/{sample}.lr_subreads.bam"),
        subreads_ref_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.lr_subreads.srt.bam"
    threads: 16
    resources:
        mem_mb = 120*1024
    shell:
        """
        pbmm2 align -j {threads} {input.consensus_fasta} {input.lr_subreads_bam} {output.subreads_ref_unsort_bam}
        samtools sort -@ {threads} -m 4G -T {wildcards.sample} -O bam -o {output.subreads_ref_bam} {output.subreads_ref_unsort_bam}
        samtools index -@ {threads} {output.subreads_ref_bam}
        """

rule personal_ref_hifi_pbmm2:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        lr_hifi_fastqs = config['lr_hifi_fastqs']
    output:
        hifi_ref_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.lr_hifi.srt.bam"
    threads: 8
    resources:
        mem_mb = 64*1024
    shell:
        """
        pbmm2 align -j {threads} -J 1 \
            {input.consensus_fasta} \
            {input.lr_hifi_fastqs} \
            {output.hifi_ref_bam} \
            --sort --preset CCS
        """

#PAR_region
rule liftover_chrX_par_region:
    input:
        chrX_par = config["par_region"],
        chain = get_chain_input
    output:
        chrX_personal_par = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.chrX_PAR.bed"
    shell:
        """
        liftOver {input.chrX_par} {input.chain} {output.chrX_personal_par} /dev/null
        """


rule phase_assembly:
    input:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        ngs_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.srt.bam",
        subreads_ref_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.lr_subreads.srt.bam",
        hifi_ref_bam = "c6_draft_assembly/sample_assembly/{sample}/{sample}.lr_hifi.srt.bam",
        sample_switch_correct_final_vcf = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.merge.phase.switch_correct.vcf.gz",
        chrX_personal_par = "c6_draft_assembly/sample_assembly/{sample}/{sample}.personal.chrX_PAR.bed"
    output:
        hap1_fa = "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap1.fasta",
        hap2_fa = "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap2.fasta"
    params:
        sex = get_sex,
        tmp_dir = "c6_draft_assembly/sample_assembly/{sample}/{sample}_tmp",
        assembly_dir = "c6_draft_assembly/sample_assembly/{sample}/assembly"
    threads: 28
    resources:
        mem_mb = 200*1024
    shell:
        """
        python3 scripts/draft_assembly/phase_assembly.version2.py \
            --ngs-bam {input.ngs_bam} \
            --subread-bam {input.subreads_ref_bam} \
            --hifi-bam {input.hifi_ref_bam} \
            -r {input.consensus_fasta} \
            -v {input.sample_switch_correct_final_vcf} \
            -p {wildcards.sample} \
            --par-bed {input.chrX_personal_par} \
            -g {params.sex} \
            -o {params.assembly_dir} \
            -T {params.tmp_dir} \
            -t {threads}
        """

rule rename_assembly:
    input:
        hap1_fa = "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap1.fasta",
        hap2_fa = "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap2.fasta",
    output:
        hap1_rename_fa = "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap1.rename.fasta",
        hap2_rename_fa = "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap2.rename.fasta",
    threads: 1
    resources:
        mem_mb = 10*1024
    shell:
        """
        awk -v sample={wildcards.sample} '{{if(substr($1,1,1)==">") print ">"sample"_hap1_"substr($1,2,length($1));else print $0}}' {input.hap1_fa} > {output.hap1_rename_fa}
        awk -v sample={wildcards.sample} '{{if(substr($1,1,1)==">") print ">"sample"_hap2_"substr($1,2,length($1));else print $0}}' {input.hap2_fa} > {output.hap2_rename_fa}
        """
