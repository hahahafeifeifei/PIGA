rule all_phase_snv:
    input:
        expand("c4_phase_snv/sample_vcf/{sample}/{sample}.shapeit.vcf.gz", sample = samples_list)

rule generate_sample_vcf_to_phase:
        input:
            merge_merfin_filter_vcf = get_merge_merfin_filter_vcf_input,
            ref = config['reference']['CHM13']
        output:
            unphase_sample_vcf = "c4_phase_snv/sample_vcf/{sample}/{sample}.consensus.vcf.gz"
        threads: 4
        shell:
            """
            bcftools view --threads {threads} -s {wildcards.sample} {input.merge_merfin_filter_vcf} | \
            bcftools view --threads {threads} -e 'GT=="0/0" | GT=="0|0"' | \
            bcftools norm --threads {threads} -m +any -f {input.ref} | \
            whatshap unphase - | bgzip -@ {threads} -c > {output.unphase_sample_vcf}
            """

rule consensus_vcf_whatshap_phase:
    input:
        ref = config['reference']['CHM13'],
        unphase_sample_vcf = "c4_phase_snv/sample_vcf/{sample}/{sample}.consensus.vcf.gz",
        zmw_bam = get_zmw_bam_input
    output:
        phase_sample_vcf = "c4_phase_snv/sample_vcf/{sample}/{sample}.consensus.whatshap.vcf.gz"
    threads: 1
    resources:
        max_mem_gb = 16
    shell:
        """
        whatshap phase --ignore-read-groups \
            --reference {input.ref} \
            --output {output.phase_sample_vcf} \
            {input.unphase_sample_vcf} \
            {input.zmw_bam}

        tabix {output.phase_sample_vcf}
        """

rule merge_consensus_whatshap_vcf:
    input:
        phase_sample_vcfs = expand("c4_phase_snv/sample_vcf/{sample}/{sample}.consensus.whatshap.vcf.gz", sample = samples_list),
        ref = config['reference']['CHM13']
    output:
        merge_whatshap_vcf = f"c4_phase_snv/merged_vcf/{config['prefix']}.consensus.whatshap.vcf.gz",
        merge_whatshap_filter_vcf = f"c4_phase_snv/merged_vcf/{config['prefix']}.consensus.whatshap.unphase_singleton_filter.vcf.gz"
    threads: 16
    shell:
        """
        bcftools merge --threads {threads} -0 -m any {input.phase_sample_vcfs} -O z -o {output.merge_whatshap_vcf}
        
        bcftools norm --threads {threads} -m -any -f {input.ref} {output.merge_whatshap_vcf} | \
        bcftools view --threads {threads} -v snps -e "N_PASS(GT=='0|1' || GT=='1|0')==0 && (AN-AC<=1 || AC<=1)" -o {output.merge_whatshap_filter_vcf}
        tabix -f {output.merge_whatshap_filter_vcf}
        """

# This step requires at least 20 samples in the input vcf file. Or shapeit4 will not run successfully.
rule chr_consensus_vcf_shapeit:
    input:
        merge_whatshap_filter_vcf = f"c4_phase_snv/merged_vcf/{config['prefix']}.consensus.whatshap.unphase_singleton_filter.vcf.gz",
        genetic_map = config['genetic_map']
    output:
        consensus_whatshap_shapeit_vcf = f"c4_phase_snv/chr_vcf/{{chr}}/{config['prefix']}.consensus.whatshap.shapeit.{{chr}}.vcf.gz"
    threads: 16
    resources:
        max_mem_gb = 100
    shell:
        """
        shapeit4 --input {input.merge_whatshap_filter_vcf} \
            --map {input.genetic_map} \
            --region {wildcards.chr} \
            --pbwt-depth 8 -T {threads} --sequencing --use-PS 0.0001 \
            --out {output.consensus_whatshap_shapeit_vcf}
        tabix -f {output.consensus_whatshap_shapeit_vcf}
        """

rule concat_consensus_vcf_shapeit:
    input:
        consensus_whatshap_shapeit_vcfs = expand("c4_phase_snv/chr_vcf/{chr}/{prefix}.consensus.whatshap.shapeit.{chr}.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"], prefix = config['prefix'])
    output:
        concat_consensus_whatshap_shapeit_vcf = f"c4_phase_snv/merged_vcf/{config['prefix']}.consensus.whatshap.shapeit.vcf.gz"
    threads: 16
    shell:
        """        
        bcftools concat {input.consensus_whatshap_shapeit_vcfs} \
            --threads {threads} -Oz -o {output.concat_consensus_whatshap_shapeit_vcf}
        """


rule prepare_sample_consensus_vcf:
    input:
        concat_consensus_whatshap_shapeit_vcf = f"c4_phase_snv/merged_vcf/{config['prefix']}.consensus.whatshap.shapeit.vcf.gz"
    output:
        sample_vcf = "c4_phase_snv/sample_vcf/{sample}/{sample}.shapeit.vcf.gz"
    threads: 4
    shell:
        """
        bcftools view --threads {threads} -s {wildcards.sample} {input.concat_consensus_whatshap_shapeit_vcf} | \
            bcftools view --threads {threads} -i "GT!='RR'" -o {output.sample_vcf}
        tabix -f {output.sample_vcf}
        """
