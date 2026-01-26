rule all_generate_personal_reference:
    input:
        expand("c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta", sample=samples_list),
        expand("c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.chain", sample=samples_list)

rule external_pangenome_index:
    input:
        external_pangenome = config['external_pangenome']
    output:
        pangenome_gbz = f"c5_personal_ref/external_pangenome/external.gbz",
        pangenome_dist = f"c5_personal_ref/external_pangenome/external.dist",
        pangenome_ri = f"c5_personal_ref/external_pangenome/external.ri",
        pangenome_hapl = f"c5_personal_ref/external_pangenome/external.hapl"
    threads: 16
    resources:
        mem_mb = 200*1024
    shell:
        """
        cp {input.external_pangenome} {output.pangenome_gbz}
        vg index -t {threads} -j {output.pangenome_dist} --no-nested-distance {output.pangenome_gbz}
        vg gbwt --num-threads {threads} -r {output.pangenome_ri} -Z {output.pangenome_gbz}
        vg haplotypes -t {threads} -H {output.pangenome_hapl} {output.pangenome_gbz}
        """
    
rule personal_pangenome:
    input:
        pangenome_gbz = f"c5_personal_ref/external_pangenome/external.gbz",
        pangenome_hapl = f"c5_personal_ref/external_pangenome/external.hapl",
        sr_fq1 = config['sr_fastqs'][0],
        sr_fq2 = config['sr_fastqs'][1],
        lr_hifi_fastqs = config['lr_hifi_fastqs']
    output:
        kmer_fq_list = "c5_personal_ref/sample_reference/{sample}/{sample}.fq.list",
        sample_kff = "c5_personal_ref/sample_reference/{sample}/{sample}.kff",
        sample_chop_gbz = "c5_personal_ref/sample_reference/{sample}/{sample}.chop.gbz",
        sample_gfa = "c5_personal_ref/sample_reference/{sample}/{sample}.gfa",
        sample_gbz = "c5_personal_ref/sample_reference/{sample}/{sample}.giraffe.gbz",
        sample_min = "c5_personal_ref/sample_reference/{sample}/{sample}.shortread.withzip.min",
        sample_dist = "c5_personal_ref/sample_reference/{sample}/{sample}.dist",
        sample_zipcodes = "c5_personal_ref/sample_reference/{sample}/{sample}.shortread.zipcodes"
    params:
        prefix = "c5_personal_ref/sample_reference/{sample}/{sample}"
    threads: 8
    resources:
        mem_mb = 200*1024
    shell:
        """
        ls {input.sr_fq1} {input.sr_fq2} {input.lr_hifi_fastqs} > {output.kmer_fq_list}
        mkdir -p {params.prefix}_tmp
        kmc -k29 -okff -t{threads} -hp @{output.kmer_fq_list} {params.prefix} {params.prefix}_tmp
        rm -rf {params.prefix}_tmp

        vg haplotypes -t {threads} --include-reference --num-haplotypes 8 --linear-structure -i {input.pangenome_hapl} -k {output.sample_kff} -g {output.sample_chop_gbz} {input.pangenome_gbz}
        vg view {output.sample_chop_gbz} | vg mod -t {threads} -u - | vg view - > {output.sample_gfa}
        vg autoindex -g {output.sample_gfa} -w giraffe -t {threads} -p {params.prefix}
        """

rule graph_call:
    input:
        lr_zmw_fastqs = config['lr_zmw_fastqs'],
        sr_fq1 = config['sr_fastqs'][0],
        sr_fq2 = config['sr_fastqs'][1],
        sample_gfa = "c5_personal_ref/sample_reference/{sample}/{sample}.gfa",
        sample_gbz = "c5_personal_ref/sample_reference/{sample}/{sample}.giraffe.gbz",
        sample_min = "c5_personal_ref/sample_reference/{sample}/{sample}.shortread.withzip.min",
        sample_dist = "c5_personal_ref/sample_reference/{sample}/{sample}.dist",
        sample_zipcodes = "c5_personal_ref/sample_reference/{sample}/{sample}.shortread.zipcodes"
    output:
        sample_gam = "c5_personal_ref/sample_reference/{sample}/{sample}.gam",
        sample_pack = "c5_personal_ref/sample_reference/{sample}/{sample}.pack",
        sample_vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.vcf.gz"
    threads: 8
    resources:
        mem_mb = 100 * 1024,
    shell:
        """
        GraphAligner -t {threads} -g {input.sample_gfa} -f {input.lr_zmw_fastqs} -a {output.sample_gam} \
        --seeds-mxm-length 30 --seeds-mem-count 10000 --bandwidth 15 --multimap-score-fraction 0.99 --min-alignment-score 100 --clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 --max-trace-count 5 --mem-index-no-wavelet-tree
        vg giraffe -t {threads} --named-coordinates -Z {input.sample_gbz} -m {input.sample_min} -z {input.sample_zipcodes} -d {input.sample_dist} -f {input.sr_fq1} -f {input.sr_fq2} >> {output.sample_gam}
        vg pack -t {threads} -x {input.sample_gfa} -g {output.sample_gam} -o {output.sample_pack}
        vg call -t {threads} {input.sample_gfa} -k {output.sample_pack} -s {wildcards.sample} -S CHM13 | bgzip -c > {output.sample_vcf}
        tabix -f {output.sample_vcf}
        """

rule gam_call_variants_filter:
    input:
        vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.vcf.gz"
    output:
        filtered_vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.filter.vcf.gz"
    threads: 2
    resources:
        mem_mb = 20*1024
    shell:
        """
        dp_mean=$(bcftools query -f "%DP\\n" {input.vcf} | awk '{{if($1>=100)print 100;else print$0}}' | csvtk -H -t summary -f 1:mean)
        dp_sd=$(bcftools query -f "%DP\\n" {input.vcf} | awk '{{if($1>=100)print 100;else print$0}}' | csvtk -H -t summary -f 1:stdev)
        bcftools view -f PASS -i "GT='AA' && FORMAT/DP<=$(awk -v dp_sd=${{dp_sd}} -v dp_mean=${{dp_mean}} 'BEGIN{{print dp_mean+3*dp_sd}}' )" {input.vcf} | \
        bcftools sort -Oz -o {output.filtered_vcf}
        tabix -f {output.filtered_vcf}
        """

rule filtered_variants_ref_consensus:
    input:
        vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.filter.vcf.gz",
        ref = config['reference']['CHM13']
    output:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        chain = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.chain"
    resources:
        mem_mb = 20*1024
    shell:
        """
        bcftools consensus -f {input.ref} -H 1 -c {output.chain} {input.vcf} > {output.consensus_fasta}
        """

