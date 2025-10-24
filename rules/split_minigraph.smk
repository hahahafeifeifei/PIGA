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
        return config["hap1_adaptor_masked_fa"][wildcards.sample]
    else:
        return "c6_draft_assembly/result/{sample}/assembly/{sample}.hap1.adaptor_masked.fasta"

def get_hap2_adaptor_masked_fa_input(wildcards):
    if "hap2_adaptor_masked_fa" in config:
        return config["hap2_adaptor_masked_fa"][wildcards.sample]
    else:
        return "c6_draft_assembly/result/{sample}/assembly/{sample}.hap2.adaptor_masked.fasta"


#TODO: use adaptor_masked fa or not?
rule split_fa_by_chr:
    input:
        hap1_adaptor_masked_fa = get_hap1_adaptor_masked_fa_input,
        hap2_adaptor_masked_fa = get_hap2_adaptor_masked_fa_input
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
        seqfile = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.seqfile"
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
        seqfile = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.seqfile"
    output:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.gfa"
    log:
        "logs/minigraph_by_chr/{chr}.log"
    threads: 16
    shell:
        """

        awk -v OFS='\t' '{{print $2}}' {input.seqfile} > c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.minigraph.seqfile
        
        minigraph -c -x ggs -l 10000 \
        -t {threads} \
        chr_gfa/t2t.grch38.58hifi.{chr}.gfa \
        $(cat c7_graph_construction/chr_mc/{chr}/{config['prefix']}.{chr}.minigraph.seqfile) \
        > {output.gfa}
        
        """

        
rule minigraph_vg_snarls:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.gfa"
    output:
        snarls = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.snarls",
        snarls_txt = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.snarls.txt"
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
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.gfa"
    output:
        og = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.og"
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
        gfa=f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.gfa",
        snarls_txt=f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.snarls.txt",
        og=f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.og"
    output:
        node_split_gfa=f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.node_split.gfa",
        node_edge_split_gfa=f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.node_split.edge_split.gfa",
        info=f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.node_split.info",
        node_split_rgfa = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.node_split.rgfa"
    log:
        "logs/border_node_select/{chr}.log"
    shell:
        """
        python3 script/graph-construct/gfa_border_node_select.py \
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
        node_split_rgfa = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.node_split.rgfa"
    output:
        fa = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.fa.gz",
        paf = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.paf",
        gaf = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.gaf.gz",
        filter_gaf = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.filter.gaf"
    params:
        mapCores = 4,
        outputGAFDir = "c7_graph_construction/chr_mc/{chr}/graphmap",
        reference = "CHM13",
        delFilter = 10000000,
        logFile = "c7_graph_construction/chr_mc/{chr}/graphmap/graphmap.log",
        tmpDir = "c7_graph_construction/chr_mc/{chr}/tmp_cactus"
    threads: 16
    shell:
        """
        rm -r c7_graph_construction/chr_mc/{wildcards.chr}/jobstore

        cactus-graphmap c7_graph_construction/chr_mc/{wildcards.chr}/jobstore {input.node_split_rgfa} {output.paf} --outputGAFDir {params.outputGAFDir} --outputFasta {output.fa} --reference {params.reference} --mapCores {params.mapCores} --delFilter {params.delFilter} --defaultPreemptable --maxNodes {threads} --logFile {params.logFile} --workDir {params.tmpDir}

        zcat {output.gaf} | \
            gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 \
            > {output.filter_gaf}
        
        """

# vg trunk will generate a list of subgraph.gfa files. Here I use subgraph_0.gfa as output to make sure the rule can work.
rule vg_chunk:
    input:
        node_edge_split_gfa = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.edge_split.gfa"
    output:
        subgraph_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/{config['prefix']}.{{chr}}.subgraph_0.gfa"
    log:
        "logs/vg_chunk/{chr}.log"
    params:
        prefix = config['prefix']
    shell:
        """
        mkdir c7_graph_construction/chr_mc/{wildcards.chr}/subgraph
        vg chunk -C -x {input.node_edge_split_gfa} --prefix c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/{params.prefix}.{wildcards.chr}.subgraph -O gfa
        
        """

rule generate_subgraph_bed:
    input:
        subgraph_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/{config['prefix']}.{{chr}}.subgraph_0.gfa",
        paf = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.paf",
        gaf = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.filter.gaf"
    output:
        subgraph_paf = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/{config['prefix']}.{{chr}}.subgraph_0.paf",
        subgraph_bed = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/{config['prefix']}.{{chr}}.subgraph_0.bed"
    log:
        "logs/generate_subgraph_bed/{chr}.log"
    params:
        prefix = config['prefix']
    shell:
        """
        python3 scripts/graph-construct/subgraph_paf_fa.py \
        c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/{params.prefix}.{wildcards.chr}.subgraph \
        {input.paf} \
        {input.gaf} \
        c7_graph_construction/chr_mc/{wildcards.chr}/subgraph
        """

#TODO: here idk why it is not ok if I replace "t2t.grch38.58hifi.1064zmw" by config['prefix']
rule process_subgraph:
    input:
        seqfile = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.seqfile",
        info = f"c7_graph_construction/chr_mc/{{chr}}/{config['prefix']}.{{chr}}.minigraph.node_split.info"
    output:
        subgraph_seqfile = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph0/CKCG.{chr}.subgraph_0.seqfile"
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
        for chr in `ls -d chr*`;do \
          for subgraph_id in `ls -d ${{chr}}/subgraph/subgraph*`;do idx=$(echo ${{subgraph_id}} | awk '{{split($1,a,"/");print a[length(a)]}}' | awk -v FS='subgraph' '{{print $2}}'); echo -e ${{chr}}"\\t"${{idx}}; done;
        done > c7_graph_construction/chr_mc/subgraph.list
        """
