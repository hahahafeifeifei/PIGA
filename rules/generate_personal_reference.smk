rule all_generate_personal_reference:
    input:
        expand("c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta", sample=samples_list),
        expand("c5_personal_ref/consensus_fasta/{sample}/{sample}.personal_ref.chain", sample=samples_list)

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
        max_mem_gb = 200
    shell:
        """
        cp {input.external_pangenome} {output.gbz}
        vg index -t {threads} -j {output.dist} --no-nested-distance {output.gbz}
        vg gbwt --num-threads {threads} -r {output.ri} -Z {output.gbz}
        vg haplotypes -t {threads} -H {output.hapl} {output.gbz}
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
        sample_gbz = "c5_personal_ref/sample_reference/{sample}/{sample}.gbz",
        sample_min = "c5_personal_ref/sample_reference/{sample}/{sample}.min",
        sample_dist = "c5_personal_ref/sample_reference/{sample}/{sample}.dist",
        sample_zipcodes = "c5_personal_ref/sample_reference/{sample}/{sample}.zipcodes"
    params:
        prefix = "c5_personal_ref/sample_reference/{sample}/{sample}"
    threads: 8
    resources:
        max_mem_gb = 150
    shell:
        """
        ls {input.sr_fq1} {input.sr_fq2} {input.lr_hifi_fastqs} > {output.kmer_fq_list}
        mkdir {params.prefix}_tmp
        kmc -k29 -m{resource.max_mem_gb} -okff -t{threads} -hp @{output.kmer_fq_list} {params.prefix} {params.prefix}_tmp
        rm -rf {params.prefix}_tmp

        vg haplotypes -t {threads} --include-reference --num-haplotypes 8 --linear-structure -i {input.pangenome_hapl} -k {output.sample_kff} -g {output.sample_chop_gbz} {input.pangenome_gbz}
        vg view {output.sample_chop_gbz} | vg mod -u - | vg view - > {output.sample_gfa}
        vg autoindex -g {output.sample_gfa} -w giraffe -t {threads} -p {params.prefix}
        """

rule graph_call:
    input:
        lr_zmw_fastqs = config['lr_zmw_fastqs'],
        sr_fq1 = config['sr_fastqs'][0],
        sr_fq2 = config['sr_fastqs'][1],
        sample_gfa = "c5_personal_ref/sample_reference/{sample}/{sample}.gfa",
        sample_gbz = "c5_personal_ref/sample_reference/{sample}/{sample}.gbz",
        sample_min = "c5_personal_ref/sample_reference/{sample}/{sample}.min",
        sample_dist = "c5_personal_ref/sample_reference/{sample}/{sample}.dist",
        sample_zipcodes = "c5_personal_ref/sample_reference/{sample}/{sample}.zipcodes"
    output:
        sample_gam = "c5_personal_ref/sample_reference/{sample}/{sample}.gam",
        sample_pack = "c5_personal_ref/sample_reference/{sample}/{sample}.pack",
        sample_vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.vcf"
    threads: 8
    resources: 
        max_mem_gb = 100
    shell:
        """
        GraphAligner -t {threads} -g {input.sample_gfa} -f {input.lr_zmw_fastqs} -a {output.sample_gam} -x vg --multimap-score-fraction 1
        vg giraffe --named-coordinates -Z {input.sample_gbz} -m {input.sample_min} -z {input.sample_zipcodes} -d {input.sample_dist} -f {input.sr_fq1} -f {input.sr_fq2} >> {output.sample_gam}
        vg pack -t {threads} -x {input.sample_gfa} -g {output.sample_gam} -o {ouput.sample_pack}
        vg call -t {threads} {input.sample_gfa} -k {ouput.sample_pack} -s {wildcards.sample} > {output.sample_vcf}
        """

rule gam_call_variants_filter:
    input:
        vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.vcf"
    output:
        filtered_vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.filter.vcf.gz"
    threads: 2
    shell:
        """
        dp_mean=$(bcftools query -f "%DP\n" {input.vcf} | awk '{{if($1>=100)print 100;else print$0}}' | csvtk -H -t summary -f 1:mean)                            
        dp_sd=$(bcftools query -f "%DP\n" {input.vcf} | awk '{{if($1>=100)print 100;else print$0}}' | csvtk -H -t summary -f 1:stdev)
        bcftools view -f PASS -i "GT='AA' && FORMAT/DP<=$(awk -v dp_sd=${{dp_sd}} -v dp_mean=${{dp_mean}} 'BEGIN{{print dp_mean+3*dp_sd}}' )" {input.vcf} | \
        bcftools sort -o {output.filtered_vcf}
        tabix -f {output.filtered_vcf}
        """

rule filtered_variants_ref_consensus:
    input:
        vcf = "c5_personal_ref/sample_reference/{sample}/{sample}.filter.vcf.gz",
        ref = config['reference']['CHM13']
    output:
        consensus_fasta = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta",
        chain = "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.chain"
    shell:
        """
        bcftools consensus -f {input.ref} -H 1 -c {output.chain} {input.vcf} > {output.consensus_fasta}
        """

