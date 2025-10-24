rule all_phase_snv:
    input:
        expand("c4_phase_snv/shapeit4/samples/{sample}/{sample}.shapeit4.vcf.gz", sample = config['samples'])


rule generate_sample_vcf_to_phase:
        input:
            merge_merfin_filter_vcf = get_merge_merfin_filter_vcf_input,
            ref = config['reference']['CHM13'],
            # sr_scaffold_vcf = get_sr_scaffold_vcf_input
        output:
            unphase_sample_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.consensus.analysis_set.vcf.gz",
            # sample_sr_scaffold_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.srs_scaffold.analysis_set.vcf.gz"
        threads: 1
        shell:
            """
            bcftools view --threads {threads} -s {wildcards.sample} {input.merge_merfin_filter_vcf}| bcftools view --threads {threads} -e 'GT=="0/0" | GT=="0|0"' | bcftools norm --threads {threads} -m +any -f {input.ref} | whatshap unphase - | bgzip -@ {threads} -c > {output.unphase_sample_vcf}

            # bcftools view --threads {threads} -s {wildcards.sample}-WGS {input.sr_scaffold_vcf} | bcftools view --threads {threads} -e 'GT=="0/0" | GT=="0|0"' | sed "s/{wildcards.sample}-WGS/{wildcards.sample}/g" | bgzip -@ {threads} -c > {output.sample_sr_scaffold_vcf}
            """

rule consensus_vcf_whatshap_phase:
    input:
        ref = config['reference']['CHM13'],
        unphase_sample_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.consensus.analysis_set.vcf.gz",
        # sample_sr_scaffold_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.srs_scaffold.analysis_set.vcf.gz",
        zmw_bam = get_zmw_bam_input
    output:
        phase_sample_vcf = "c4_phase_snv/whatshap/{sample}/{sample}.consensus.whatshap.analysis_set.vcf.gz"
    threads: 1
    resources:
        #16G
        mem_mb = 16000
    shell:
        """
      #   whatshap phase --ignore-read-groups \
      #       --reference {input.ref} \
      #       --output {output.phase_sample_vcf} \
      #       {input.unphase_sample_vcf} \
      #       {input.zmw_bam} {input.sample_sr_scaffold_vcf} 
        
        whatshap phase --ignore-read-groups \
            --reference {input.ref} \
            --output {output.phase_sample_vcf} \
            {input.unphase_sample_vcf} \
            {input.zmw_bam}

        tabix {output.phase_sample_vcf}
        """

rule merge_consensus_whatshap_vcf:
    input:
        phase_sample_vcfs = expand("c4_phase_snv/whatshap/{sample}/{sample}.consensus.whatshap.analysis_set.vcf.gz", sample = config['samples']),
        ref = config['reference']['CHM13']
    output:
        phase_sample_vcf_list = "c4_phase_snv/whatshap/merge/vcf.list",
        merge_whatshap_vcf = f"c4_phase_snv/whatshap/merge/{config['prefix']}.consensus.whatshap.analysis_set.vcf.gz",
        merge_whatshap_filter_vcf = f"c4_phase_snv/whatshap/merge/{config['prefix']}.consensus.whatshap.unphase_singleton_filter.analysis_set.vcf.gz"
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
        merge_whatshap_filter_vcf = f"c4_phase_snv/whatshap/merge/{config['prefix']}.consensus.whatshap.unphase_singleton_filter.analysis_set.vcf.gz",
        # sr_scaffold_vcf = f"c3_merge_snv/callset/srs_lrs_compare/{config['prefix']}.analysis_set.call_set.consensus.srs_scaffold.vcf.gz",
        genetic_map = config['genetic_map'],
        # topmed_east_asian = "/storage/yangjianLab/wangyifei/resource/TOPMed/east_asian/merge_vcf/{chr}/TOPMed_WGS_freeze.8.east_asian.merge.{chr}.filter.snps.liftover_chm13.shapeit4.vcf.gz",
    output:
      # consensus_whatshap_shapeit4_scaffold = f"c4_phase_snv/shapeit4/{{chr}}/{config['prefix']}.consensus.whatshap.topmed_eas.scaffold.shapeit4.analysis_set.{{chr}}.vcf.gz",
        consensus_whatshap_shapeit4_vcf = f"c4_phase_snv/shapeit4/{{chr}}/{config['prefix']}.consensus.whatshap.topmed_eas.shapeit4.analysis_set.{{chr}}.vcf.gz"
    threads: 16
    resources:
        #100G
        mem_mb = 100000
    shell:
        """
        # shapeit4 --input {input.merge_whatshap_filter_vcf} \
        #     --map {input.genetic_map} \
        #     --region {wildcards.chr} \
        #     --pbwt-depth 8 -T {threads} --sequencing --use-PS 0.0001 \
        #     --reference {input.topmed_east_asian} \
        #     --out {output.consensus_whatshap_shapeit4_scaffold}
        #    --scaffold {input.sr_scaffold_vcf}

        # tabix -f {output.consensus_whatshap_shapeit4_scaffold}

        shapeit4 --input {input.merge_whatshap_filter_vcf} \
            --map {input.genetic_map} \
            --region {wildcards.chr} \
            --pbwt-depth 8 -T {threads} --sequencing --use-PS 0.0001 \
            --out {output.consensus_whatshap_shapeit4_vcf}
           # --scaffold {output.consensus_whatshap_shapeit4_scaffold}

        tabix -f {output.consensus_whatshap_shapeit4_vcf}
        """

rule concat_consensus_vcf_shapeit4:
    input:
        consensus_whatshap_shapeit4_vcfs = expand("c4_phase_snv/shapeit4/{chr}/{prefix}.consensus.whatshap.topmed_eas.shapeit4.analysis_set.{chr}.vcf.gz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX"], prefix = config['prefix'])
    output:
        consensus_whatshap_shapeit4_vcf_list = "c4_phase_snv/shapeit4/vcf.list",
        concat_consensus_whatshap_shapeit4_vcf = f"c4_phase_snv/shapeit4/{config['prefix']}.consensus.whatshap.shapeit4.analysis_set.vcf.gz"
    threads: 16
    shell:
        """
        ls c4_phase_snv/shapeit4/*/{config['prefix']}.consensus.whatshap.shapeit4.analysis_set.chr*.vcf.gz > {output.consensus_whatshap_shapeit4_vcf_list}
        
        bcftools concat -f {output.consensus_whatshap_shapeit4_vcf_list} \
            --threads {threads} \
            -Oz \
            -o {output.concat_consensus_whatshap_shapeit4_vcf}
        """


rule prepare_sample_consensus_vcf:
    input:
        concat_consensus_whatshap_shapeit4_vcf = f"c4_phase_snv/shapeit4/{config['prefix']}.consensus.whatshap.shapeit4.analysis_set.vcf.gz"
    output:
        sample_vcf = temp("c4_phase_snv/shapeit4/samples/{sample}/{sample}.shapeit4.vcf.gz")
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





