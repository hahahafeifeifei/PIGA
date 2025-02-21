# merfin polish

def get_sr_fastqs(wildcards):
    return config["sr_fastqs"][wildcards.sample]

def get_sr_unmapped_fastqs(wildcards):
    return config["sr_unmapped_fastqs"][wildcards.sample]
    
def get_Q20_pbmm2_map_input_fastqs(wildcards):
    return config["lr_Q20_fastqs"][wildcards.sample]
#use merylIndex to represent the whole {sample}.meryl directory.
rule prepare_sample_kmer:
    input:
        ngs_R1_fastq = get_sr_fastqs[0],
        ngs_R2_fastq = get_sr_fastqs[1],
        ngs_R1_unmapped_fastq = get_sr_unmapped_fastqs[0],
        ngs_R2_unmapped_fastq = get_sr_unmapped_fastqs[1],
        hifi_fastq = get_Q20_pbmm2_map_input_fastqs
    output:
        sample_meryl = "c4_phase_snv/meryl/{sample}/{sample}.meryl/merylIndex"
    threads: 4
    resources:
        #30G
        mem_mb = 30000
    params:
        mem = 30
    shell:
        """
        meryl count k=21 memory={params.mem} threads={threads} output c4_phase_snv/meryl/{wildcards.sample}/R1.meryl {input.ngs_R1_fastq}
        meryl count k=21 memory={params.mem} threads={threads} output c4_phase_snv/meryl/{wildcards.sample}/R2.meryl {input.ngs_R2_fastq}
        meryl count k=21 memory={params.mem} threads={threads} output c4_phase_snv/meryl/{wildcards.sample}/unpaired.R1.meryl {input.ngs_R1_unmapped_fastq}
        meryl count k=21 memory={params.mem} threads={threads} output c4_phase_snv/meryl/{wildcards.sample}/unpaired.R2.meryl {input.ngs_R2_unmapped_fastq}
        meryl union-sum output c4_phase_snv/meryl/{wildcards.sample}/{wildcards.sample}-WGS.meryl c4_phase_snv/meryl/{wildcards.sample}/R1.meryl c4_phase_snv/meryl/{wildcards.sample}/R2.meryl c4_phase_snv/meryl/{wildcards.sample}/unpaired.R1.meryl c4_phase_snv/meryl/{wildcards.sample}/unpaired.R2.meryl
        rm -rf c4_phase_snv/meryl/{wildcards.sample}/R1.meryl c4_phase_snv/meryl/{wildcards.sample}/R2.meryl c4_phase_snv/meryl/{wildcards.sample}/unpaired.R1.meryl c4_phase_snv/meryl/{wildcards.sample}/unpaired.R2.meryl

        meryl count k=21 memory={params.mem} threads={threads} output c4_phase_snv/meryl/{wildcards.sample}/hifi.meryl {input.hifi_fastq}
        meryl union-sum output c4_phase_snv/meryl/{wildcards.sample}/{wildcards.sample}.meryl c4_phase_snv/meryl/{wildcards.sample}/{wildcards.sample}-WGS.meryl c4_phase_snv/meryl/{wildcards.sample}/hifi.meryl

        rm -rf c4_phase_snv/meryl/{wildcards.sample}/hifi.meryl
        """

rule merfin_filter:
    input:
        consensus_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.vcf.gz",
        sample_consensus_vcf = "c3_merge_snv/samples/{sample}/{sample}.consensus.vcf",
        ref = config['CHM13'],
        ref_mer = directory("/storage/yangjianLab/wangyifei/resource/Reference/CHM13/CHM13_21mer_meryl"),
        sample_WGS_meryl = "c4_phase_snv/meryl/{sample}/{sample}-WGS.meryl/merylIndex",
        sample_consensus_miss_vcf = "c3_merge_snv/samples/{sample}/{sample}.consensus.miss.vcf.gz"
    output:
        sample_consensus_filter_vcf = "c4_phase_snv/merfin/{sample}/{sample}.consensus.filter.vcf.gz",
        sample_consensus_final_vcf = "c4_phase_snv/merfin/{sample}/{sample}.consensus.final.vcf.gz"
    threads: 6
    resources:
        #60G
        mem_mb = 60000
    params:
        mem = 60
    shell:
        """
        merfin -filter \
            -sequence {input.ref} \
            -readmers c4_phase_snv/meryl/{wildcards.sample}/{wildcards.sample}-WGS.meryl \
            -seqmers {input.ref_mer} \
            -vcf {input.consensus_vcf} \
            -output c4_phase_snv/merfin/{wildcards.sample}/{wildcards.sample}.consensus \
            -threads {threads} \
            -memory {params.mem}

        bgzip -@ {threads} c4_phase_snv/merfin/{wildcards.sample}/{wildcards.sample}.consensus.filter.vcf
        tabix -f {output.sample_consensus_filter_vcf}

        bcftools concat --threads {threads} -a {output.sample_consensus_filter_vcf} {input.sample_consensus_miss_vcf} -o {output.sample_consensus_final_vcf}
        """

rule merge_merfin_vcf:
    input:
        expand("c4_phase_snv/merfin/{sample}/{sample}.consensus.final.vcf.gz", sample = config['samples'])
    output:
        merfin_vcf_list = "c4_phase_snv/merfin/merge/vcf.list",
        merge_merfin_vcf = "c4_phase_snv/merfin/merge/CKCG.CHM13.consensus.phase1.call_set.vcf.gz",
        merge_merfin_filter_vcf = "c4_phase_snv/merfin/merge/CKCG.CHM13.consensus.phase1.call_set.hwe_missing_filter.vcf.gz",
    threads: 28
    resources:
        #50G
        mem_mb = 50000
    shell:
        """
        ls c4_phase_snv/merfin/*/*.consensus.final.vcf.gz > {output.merfin_vcf_list}

        bcftools merge -0 -m none --threads {threads} \
            -l {output.merfin_vcf_list} | \
            bcftools sort -m 50G | \
            bcftools plugin fill-tags --threads {threads} | \
            bcftools view --threads {threads} -i "AC!=0" -o {output.merge_merfin_vcf}

        bcftools view --threads {threads} -i "HWE>=1e-6 && NS>1007" {output.merge_merfin_vcf} -o {output.merge_merfin_filter_vcf}
        """

rule generate_sample_vcf_to_phase:
        input:
            merge_merfin_filter_vcf = "c4_phase_snv/merfin/merge/CKCG.CHM13.consensus.phase1.call_set.hwe_missing_filter.vcf.gz",
            ref = config['CHM13'],
            sr_scaffold_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.srs_scaffold.vcf.gz"
        output:
            unphase_sample_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.consensus.analysis_set.vcf.gz",
            sample_sr_scaffold_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.srs_scaffold.analysis_set.vcf.gz"
        shell:
            """
            bcftools view --threads {threads} -s {wildcards.sample} {input.merge_merfin_filter_vcf}| bcftools view --threads {threads} -e 'GT=="0/0" | GT=="0|0"' | bcftools norm --threads {threads} -m +any -f {input.ref} | whatshap unphase - | bgzip -@ {threads} -c > {output.unphase_sample_vcf}

            bcftools view --threads {threads} -s {wildcards.sample}-WGS {input.sr_scaffold_vcf} | bcftools view --threads {threads} -e 'GT=="0/0" | GT=="0|0"' | sed "s/{wildcards.sample}-WGS/{wildcards.sample}/g" | bgzip -@ {threads} -c > {output.sample_sr_scaffold_vcf}
            """

rule consensus_vcf_whatshap_phase:
    input:
        ref = config['CHM13'],
        unphase_sample_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.consensus.analysis_set.vcf.gz",
        sample_sr_scaffold_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.srs_scaffold.analysis_set.vcf.gz",
        zmw_bam = "c2_call_lr_snv/lr_mapping/{sample}/{sample}.zmw.pbmm2.bam"
    output:
        phase_sample_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.consensus.whatshap.analysis_set.vcf.gz"
    threads: 1
    resources:
        #16G
        mem_mb = 16000
    shell:
        """
        whatshap phase --ignore-read-groups \
            --reference {input.ref} \
            --output {output.phase_sample_vcf} \
            {input.unphase_sample_vcf} \
            {input.zmw_bam} {input.sample_sr_scaffold_vcf} 
        
        tabix {output.phase_sample_vcf}
        """

rule merge_consensus_whatshap_vcf:
    input:
        phase_sample_vcfs = expand("c4_phase_snv/whatshap/{sample}/{sample}.consensus.whatshap.analysis_set.vcf.gz", sample = config['samples']),
        ref = config['CHM13']
    output:
        phase_sample_vcf_list = "c4_phase_snv/whatshap/merge/vcf.list",
        merge_whatshap_vcf = "c4_phase_snv/whatshap/merge/CKCG.consensus.whatshap.analysis_set.vcf.gz",
        merge_whatshap_filter_vcf = "c4_phase_snv/whatshap/merge/CKCG.consensus.whatshap.unphase_singleton_filter.analysis_set.vcf.gz"
    threads: 16
    shell:
        """
        ls c4_phase_snv/whatshap/*/*.consensus.whatshap.analysis_set.vcf.gz > c4_phase_snv/whatshap/merge/vcf.list
        bcftools merge --threads {threads} -0 -m any -l {output.phase_sample_vcf_list} -O z -o {output.merge_whatshap_vcf}
        
        bcftools norm --threads {threads} -m -any -f {input.ref} {output.merge_whatshap_vcf} | bcftools view --threads {threads} -v snps -e "N_PASS(GT=='0|1' || GT=='1|0')==0 && (AN-AC<=1 || AC<=1)" -o {output.merge_whatshap_filter_vcf}
        tabix -f {output.merge_whatshap_filter_vcf}
        """

rule chr_consensus_vcf_shapeit4:
    input:
        merge_whatshap_filter_vcf = "c4_phase_snv/whatshap/merge/CKCG.consensus.whatshap.unphase_singleton_filter.analysis_set.vcf.gz",
        sr_scaffold_vcf = "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.srs_scaffold.vcf.gz",
        genetic_map = "/storage/yangjianLab/wangyifei/resource/genetic_map/chm13_liftover/shapeit4/genetic_map_chm13_{chr}.reform.txt",
        topmed_east_asian = "/storage/yangjianLab/wangyifei/resource/TOPMed/east_asian/merge_vcf/{chr}/TOPMed_WGS_freeze.8.east_asian.merge.{chr}.filter.snps.liftover_chm13.shapeit4.vcf.gz",
    output:
        consensus_whatshap_shapeit4_scaffold = "c4_phase_snv/shapeit4/{chr}/CKCG.consensus.whatshap.topmed_eas.scaffold.shapeit4.analysis_set.{chr}.vcf.gz",
        consensus_whatshap_shapeit4_vcf = "c4_phase_snv/shapeit4/{chr}/CKCG.consensus.whatshap.topmed_eas.shapeit4.analysis_set.{chr}.vcf.gz"
    threads: 16
    resources:
        #100G
        mem_mb = 100000
    shell:
        """
        shapeit4 --input {input.merge_whatshap_filter_vcf} \
            --map {input.genetic_map} \
            --region {wildcards.chr} \
            --pbwt-depth 8 -T {threads} --sequencing --use-PS 0.0001 \
            --scaffold {input.sr_scaffold_vcf} \
            --reference {input.topmed_east_asian} \
            --out {output.consensus_whatshap_shapeit4_scaffold}

        tabix -f {output.consensus_whatshap_shapeit4_scaffold}

        shapeit4 --input {input.merge_whatshap_filter_vcf} \
            --map {input.genetic_map} \
            --region {wildcards.chr} \
            --pbwt-depth 8 -T {threads} --sequencing --use-PS 0.0001 \
            --scaffold {output.consensus_whatshap_shapeit4_scaffold} \
            --out {output.consensus_whatshap_shapeit4_vcf}

        tabix -f {output.consensus_whatshap_shapeit4_vcf}
        """

rule concat_consensus_vcf_shapeit4:
    input:
        consensus_whatshap_shapeit4_vcfs = expand("c4_phase_snv/shapeit4/{chr}/CKCG.consensus.whatshap.topmed_eas.shapeit4.analysis_set.{chr}.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"])
    output:
        consensus_whatshap_shapeit4_vcf_list = "c4_phase_snv/shapeit4/vcf.list",
        concat_consensus_whatshap_shapeit4_vcf = "c4_phase_snv/shapeit4/CKCG.consensus.whatshap.shapeit4.analysis_set.vcf.gz"
    threads: 16
    shell:
        """
        ls c4_phase_snv/shapeit4/*/CKCG.consensus.whatshap.shapeit4.analysis_set.chr*.vcf.gz > {output.consensus_whatshap_shapeit4_vcf_list}
        
        bcftools concat -f {output.consensus_whatshap_shapeit4_vcf_list} \
            --threads {threads} \
            -Oz \
            -o {output.concat_consensus_whatshap_shapeit4_vcf}
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





