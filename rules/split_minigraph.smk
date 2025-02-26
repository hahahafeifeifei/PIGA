# chr_list= [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

def get_all_subgraph_seqfiles(wildcards, prefix):
    subgraph_combination = checkpoints.prepare_subgraph_list.get().output[0]
    pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, subgraph_id = line.strip().split("\t")
            subgraph_seqfile = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqfile"
            pairs.append(subgraph_seqfile)
    return pairs

rule all_split_minigraph:
    input:
        partial(get_all_subgraph_seqfiles, prefix = config['prefix'])
        

def get_hap1_adaptor_masked_fa_input(wildcards):
    if "hap1_adaptor_masked_fa" in config:
        return config["hap1_adaptor_masked_fa"]
    else:
        return "c6_draft_assembly/{sample}/assembly/{sample}.hap1.adaptor_masked.fasta"

def get_hap2_adaptor_masked_fa_input(wildcards):
    if "hap2_adaptor_masked_fa" in config:
        return config["hap2_adaptor_masked_fa"]
    else:
        return "c6_draft_assembly/{sample}/assembly/{sample}.hap2.adaptor_masked.fasta"


rule split_fa_by_chr:
    input:
        hap1_adaptor_masked_fa = "c6_draft_assembly/{sample}/assembly/{sample}.hap1.adaptor_masked.fasta",
        hap2_adaptor_masked_fa = "c6_draft_assembly/{sample}/assembly/{sample}.hap2.adaptor_masked.fasta"
    output:
        chr_hap1_fasta = "c7_graph_construction/chr_fa/{chr}/{sample}.hap1.fasta",
        chr_hap2_fasta = "c7_graph_construction/chr_fa/{chr}/{sample}.hap2.fasta"
    log:
        "logs/split_fa_by_chr/{chr}/{sample}.log"
    shell:
        """

        seqkit grep -w 0 -r -p _{wildcards.chr}- {input.hap1_adaptor_masked_fa}  > {output.chr_hap1_fasta}
        seqkit grep -w 0 -r -p _{wildcards.chr}- {input.hap2_adaptor_masked_fa}  > {output.chr_hap2_fasta}
        
        
        """

#TODO: reference 出现两遍
rule generate_seqfile_by_chr:
    input:
        samples_hap1_fa = expand("c7_graph_construction/chr_fa/{chr}/{sample}.hap1.fasta", sample = config['samples'], allow_missing=True),
        samples_hap2_fa = expand("c7_graph_construction/chr_fa/{chr}/{sample}.hap2.fasta", sample = config['samples'], allow_missing=True),
        chm13_ref = config['reference']['CHM13'],
        grch38_ref = config['reference']['GRCh38']
    output:
        seqfile = "c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.seqfile"
    log:
        "logs/generate_seqfile_by_chr/{chr}.log"
    shell:
        """
        seqkit grep -w 0 -r -p "\b{wildcards.chr}\b" {input.chm13_ref} > c7_graph_construction/chr_fa/{wildcards.chr}/CHM13.{wildcards.chr}.fasta
        echo -e "CHM13\tc7_graph_construction/chr_fa/{wildcards.chr}/CHM13.{wildcards.chr}.fasta" >> {output.seqfile}
        
        seqkit grep -w 0 -r -p "\b{wildcards.chr}\b" {input.grch38_ref} > c7_graph_construction/chr_fa/{wildcards.chr}/GRCh38.{wildcards.chr}.fasta
        echo -e "GRCh38\tc7_graph_construction/chr_fa/{wildcards.chr}/GRCh38.{wildcards.chr}.fasta" >> {output.seqfile}
        
        for fasta in `ls c7_graph_construction/chr_fa/{wildcards.chr}/*.fasta`;do \
            sample=$(echo $fasta | awk '{{split($1,a,"/");print a[length(a)]}}' | cut -d '.' -f 1)
            hap=$(echo $fasta | awk '{{split($1,a,"/");print a[length(a)]}}' | cut -d '.' -f 2 | cut -d "p" -f 2)
            echo -e ${{sample}}.${{hap}}"\t"${{fasta}} >> {output.seqfile}
            done
        
        """

        
#first sort by depth. Then run the minigraph
rule minigraph_by_chr:
    input:
        seqfile = "c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.seqfile",
        depth_file = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/CKCG.zmw.depth"
    output:
        gfa = "c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.minigraph.gfa"
    log:
        "logs/minigraph_by_chr/{chr}.log"
    threads: 16
    shell:
        """
        
        awk -v OFS='\\t' '{{split($1,a,".");print a[1],$2}}' {input.seqfile} | \
        csvtk -H -t join -f 1 - {input.depth_file} | \
        sort -k 3nr | \
        awk '{{print$2}}' > c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.minigraph.seqfile
        
        minigraph -c -x ggs -l 10000 \
        -t {threads} \
        chr_gfa/t2t.grch38.58hifi.{chr}.gfa \
        $(cat c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.minigraph.seqfile) \
        > {output.gfa}
        
        """
        

        
        
rule extract_minigraph_call_region:
    input:
        gfa = "c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.minigraph.gfa"
    output:
        call_region_bed = "c7_graph_construction/chr_mc/{chr}/minigraph_cov/minigraph_call.merge.bed"
    log:
        "logs/extract_minigraph_call_region/{chr}.log"
    threads: 16
    shell:
        """
        for fasta in chr_fa/{chr}/*fasta;
        do \
            sample=$(echo $fasta | cut -d "/" -f 3 | cut -d "." -f 1-2)
            minigraph -c -x asm \
            --call --vc \
            -t {threads} {input.gfa} ${{fasta}} | python3 ~/software/script/graph-construct/minigraph_call_bed_transform.py - ${{sample}} > c7_graph_construction/chr_mc/{chr}/minigraph_cov/${sample}.minigraph_call.bed
        done
        
        cat c7_graph_construction/chr_mc/{chr}/minigraph_cov/*.minigraph_call.bed > {output.call_region_bed}
        """
        
rule minigraph_cov_filter:
    input:
        gfa = "c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.minigraph.gfa",
        call_region_bed = "c7_graph_construction/chr_mc/{chr}/minigraph_cov/minigraph_call.merge.bed"
    output:
        filtered_gfa = "c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.filter.rmtips.unchop.gfa"
    log:
        "logs/minigraph_cov_filter/{chr}.log"
    shell:
        """
        python3 ~/software/script/graph-construct/minigraph_cov_filter.py \
        {input.call_region_bed} \
        {input.gfa} | \
        gfatools view -R {chr} - > c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.filter.gfa
        
        python3 ~/software/script/graph-construct/remove_gfa_tips.py \
        c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.filter.gfa \
        > c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.filter.rmtips.gfa
        
        vg convert -g c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.filter.rmtips.gfa | vg mod -u - | vg view - > {output.filtered_gfa}
        """
        
rule minigraph_vg_snarls:
    input:
        gfa = "c7_graph_construction/chr_mc/{config['prefix']}.{chr}.filter.rmtips.unchop.gfa"
    output:
        snarls = "c7_graph_construction/chr_mc/{config['prefix']}.{chr}.filter.rmtips.unchop.snarls",
        snarls_txt = "c7_graph_construction/chr_mc/{config['prefix']}.{chr}.filter.rmtips.unchop.snarls.txt"
    log:
        "logs/minigraph_vg_snarls/{chr}.log"
    shell:
        """
        vg snarls {input.gfa} > {output.snarls}
        vg view -R {output.snarls} > {output.snarls_txt}
        """

#TODO:check if the tools can run.(odgi)
rule odgi_build:
    input:
        gfa = "c7_graph_construction/chr_mc/{config['prefix']}.{chr}.filter.rmtips.unchop.gfa"
    output:
        og = "c7_graph_construction/chr_mc/{config['prefix']}.{chr}.filter.rmtips.unchop.og"
    log:
        "logs/odgi_build/{chr}.log"
    threads: 16
    shell:
        """
        source ~/miniconda3/bin/activate graph-tools
        export LD_PRELOAD=~/miniconda3/envs/graph-tools/lib/libjemalloc.so.2
        odgi build -t {threads} \
        -g {input.gfa} \
        -o {output.og}
        """
        
rule border_node_select:
    input:
        gfa="c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.gfa",
        snarls_txt="c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.snarls.txt",
        og="c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.og"
    output:
        node_split_gfa="c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.gfa",
        node_edge_split_gfa="c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.edge_split.gfa",
        info="c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.info",
        node_split_rgfa = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.rgfa"
    log:
        "logs/border_node_select/{chr}.log"
    shell:
        """
        python3 ~/software/script/graph-construct/border_node_select.py \
        {input.gfa} \
        {input.snarls_txt} \
        {input.og} \
        {output.node_split_gfa} \
        {output.node_edge_split_gfa} \
        {output.info}
        
        vg convert -g {output.node_split_gfa} -f -Q chr | \
        awk '{{if(substr($0,1,1)=="S" && $6!="SR:i:0")print$0"\tSN:Z:Other\tSO:i:0\tSR:i:1";else print$0}}' | \
        awk -v OFS='\\t' '{{if(substr($0,1,1)=="S")print$1,"s"$2,$3,$4,$5,$6;else {{if(substr($0,1,1)=="L")print$1,"s"$2,$3,"s"$4,$5,$6;else print$0}} }}' \
        > {output.node_split_rgfa}
        """
        
rule cactus_graphmap:
    input:
        node_split_rgfa = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.rgfa"
    output:
        fa = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.fa.gz",
        paf = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.paf",
        gaf = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.gaf.gz"
    params:
        mapCores = 4,
        outputGAFDir = "c7_graph_construction/chr_mc/{chr}/01.graphmap",
        reference = "CHM13",
        delFilter = 10000000,
        logFile = "c7_graph_construction/chr_mc/{chr}/01.graphmap/01.graphmap.log",
        tmpDir = "c7_graph_construction/chr_mc/{chr}/tmp_cactus"
    threads: 16
    shell:
        """
        rm -r c7_graph_construction/chr_mc/{wildcards.chr}/jobstore

        cactus-graphmap c7_graph_construction/chr_mc/{wildcards.chr}/jobstore {input.node_split_rgfa} {output.paf} --outputGAFDir {params.outputGAFDir} --outputFasta {output.fa} --reference {params.reference} --mapCores {params.mapCores} --delFilter {params.delFilter} --defaultPreemptable --maxNodes {threads} --logFile {params.logFile} --workDir {params.tmpDir}
        """
        
rule recover_29mb_edge:
    input:
        info = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.info",
        node_edge_split_gfa = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.edge_split.gfa"
    output:
        recover20mb_gfa = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.edge_split.recover20mb.gfa"
    log:
        "logs/recover_29mb_edge/{chr}.log"
    shell:
        """
        awk '{{if(NR%4!=1)print $5; if(NR%4!=0)printf $9"\\t"}}' {input.info} | \
        awk -v OFS='\\t' '{{if(NF!=1)print "L",$1,"+",$2,"+","0M"}}' | \
        cat {input.node_edge_split_gfa} - > {output.recover20mb_gfa}
        """

# vg trunk will generate a list of subgraph.gfa files. Here I use subgraph_0.gfa as output to make sure the rule can work.
rule vg_chunk:
    input:
        recover20mb_gfa = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.edge_split.recover20mb.gfa"
    output:
        subgraph_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/{config['prefix']}.{chr}.subgraph_0.gfa"
    log:
        "logs/vg_chunk/{chr}.log"
    params:
        prefix = config['prefix']
    shell:
        """
        mkdir c7_graph_construction/chr_mc/{wildcards.chr}/subgraph
        vg chunk -C -x {input.recover20mb_gfa} --prefix c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/{params.prefix}.{wildcards.chr}.subgraph -O gfa
        
        """

rule generate_subgraph_bed:
    input:
        subgraph_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/{config['prefix']}.{chr}.subgraph_0.gfa",
        paf = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.paf",
        gaf = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.gaf"
    output:
        subgraph_paf = "c7_graph_construction/chr_mc/{chr}/subgraph/{config['prefix']}.{chr}.subgraph_0.paf",
        subgraph_bed = "c7_graph_construction/chr_mc/{chr}/subgraph/{config['prefix']}.{chr}.subgraph_0.bed"
    log:
        "logs/generate_subgraph_bed/{chr}.log"
    params:
        prefix = config['prefix']
    shell:
        """
        python3 scripts/graph-construct/subgraph_paf_fa.py \
        c7_graph_construction/chr_mc/{chr}/subgraph/{params.prefix}.{chr}.subgraph \
        {input.paf} \
        {input.gaf} \
        c7_graph_construction/chr_mc/{chr}/subgraph
        """

#TODO: here idk why it is not ok if I replace "t2t.grch38.58hifi.1064zmw" by config['prefix']
rule process_subgraph:
    input:
        seqfile = "c7_graph_construction/chr_mc/{chr}/t2t.grch38.58hifi.1064zmw.{chr}.seqfile",
        info = "c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.filter.rmtips.unchop.node_split.info"
    output:
        subgraph_seqfile = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph0/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_0.seqfile"
    log:
        "logs/process_subgraph/{chr}.log"
    params:
        prefix = config['prefix'],
        chr_dir = "c7_graph_construction/chr_mc/{chr}"
    shell:
        """
        row=$[$(cat {input.info} | wc -l)-1]
        for i in `seq 0 $row`;
        do
            mkdir {params.chr_dir}/subgraph/subgraph${{i}}
            mv {params.chr_dir}/subgraph/{params.prefix}.{wildcards.chr}.subgraph_${{i}}.* {params.chr_dir}/subgraph/subgraph${{i}}
            
            mkdir {params.chr_dir}/subgraph/subgraph${{i}}/fa
            > {params.chr_dir}/subgraph/subgraph${{i}}/{params.prefix}.{wildcards.chr}.subgraph_${{i}}.seqfile

            awk 'NR>1 && NR<2509{{print$0}}' {input.seqfile} | while read line; do

                sample=$(echo $line | awk '{{print$1}}')
                fa=$(echo $line | awk '{{print$2}}')
                sample_prefix=$(echo $fa | awk '{{split($1,a,"/");split(a[length(a)],b,".fasta");print b[1]}}')
                awk -v sample=${{sample}} '{{if($1==sample)print$2}}' {params.chr_dir}/subgraph/subgraph${{i}}/{params.prefix}.{wildcards.chr}.subgraph_${{i}}.bed > {params.chr_dir}/subgraph/subgraph${{i}}/fa/${{sample_prefix}}.subgraph${{i}}.bed

                samtools faidx -r {params.chr_dir}/subgraph/subgraph${{i}}/fa/${{sample_prefix}}.subgraph${{i}}.bed ${{fa}} | awk -v sample=${{sample}} '{{if(substr($0,1,1)==">")print ">id="${{sample}}"|"substr($0,2,length($0));else print$0}}' > {params.chr_dir}/subgraph/subgraph${{i}}/fa/${{sample_prefix}}.subgraph${{i}}.fasta

                echo -e ${{sample}}"\\t"{params.chr_dir}/subgraph/subgraph${{i}}/fa/${{sample_prefix}}.subgraph${{i}}.fasta >> {params.chr_dir}/subgraph/subgraph${{i}}/{params.prefix}.{wildcards.chr}.subgraph_${{i}}.seqfile

            done
            
            awk '{{if($1=="S")print">id=_MINIGRAPH_|s"$2"\\n"$3}}' {params.chr_dir}/subgraph/subgraph${{i}}/{params.prefix}.{wildcards.chr}.subgraph_${i}.gfa > {params.chr_dir}/subgraph/subgraph${{i}}/fa/_MINIGRAPH_.subgraph${{i}}.fasta

            echo -e _MINIGRAPH_"\\t"{params.chr_dir}/subgraph/subgraph${{i}}/fa/_MINIGRAPH_.subgraph${{i}}.fasta >> {params.chr_dir}/subgraph/subgraph${{i}}/{params.prefix}.{wildcards.chr}.subgraph_${{i}}.seqfile

        done
        """


#Here subgraph_seqfile only to show the graph splitting has been finished.
checkpoint prepare_subgraph_list:
    input:
        subgraph_seqfile = expand("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph0/{prefix}.{chr}.subgraph_0.seqfile", chr=chr_list, prefix=config['prefix'])
    output:
        subgraph_list = "c7_graph_construction/chr_mc/subgraph.list" 
    shell:
        """
        for chr in `ls -d chr*`;do  for subgraph_id in `ls -d ${{chr}}/subgraph/subgraph*`;do idx=$(echo ${{subgraph_id}} | awk '{{split($1,a,"/");print a[length(a)]}}' | awk -v FS='subgraph' '{{print $2}}'); echo -e ${{chr}}"\\t"${{idx}}; done; done > c7_graph_construction/chr_mc/subgraph.list
        """
