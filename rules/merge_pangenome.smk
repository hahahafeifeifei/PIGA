
# chr_list= [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

rule all_merge_pangenome:
    input:
        f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.hapl",
        f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbz",
        expand("c7_graph_construction/graph_merge/{prefix}.{chr}.assembly.gfa", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"], prefix = config['prefix']),
        expand("c7_graph_construction/graph_merge/{prefix}.{chr}.variant.path", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"], prefix = config['prefix'])


        
###TODO:how to get the ordered subgraph list?
rule subgraph_gfa_merge:
    input:
        partial(get_all_subgraph_assembly_gfa_files, prefix = config['prefix']),
        subgraph_order_list = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/subgraph_order.list"
    output:
        merged_gfa = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gfa"
    threads: 8
    params:
        prefix = lambda wildcards: config['prefix']
    shell:
        """
        awk -v chr={wildcards.chr} '{{if($1==chr) print$2}}' {input.subgraph_order_list} | while read subgraph_id
        do
        echo c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph${{subgraph_id}}/{params.prefix}.{wildcards.chr}.subgraph_${{subgraph}}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.assembly.gfa
        done > c7_graph_construction/graph_merge/{wildcards.chr}.gfa.list
        
        python3 scripts/graph-simplification/gfa_merge.py \
        c7_graph_construction/graph_merge/{wildcards.chr}.gfa.list \
        {output.merged_gfa} \
        c7_graph_construction/graph_merge/tmp.merge.{wildcards.chr} \
        {threads}
        """

rule subgraph_variant_path_merge:
    input:
        partial(get_all_subgraph_variant_path_files, prefix = config['prefix']),
        subgraph_order_list = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/subgraph_order.list"
    output:
        merged_variant_path = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.variant.path"
    params:
        prefix = lambda wildcards: config['prefix']
    shell:
        """
        > {output.merged_variant_path}
        awk -v chr={wildcards.chr} '{{if($1==chr) print$2}}' {input.subgraph_order_list} | while read subgraph_id
        do
        cat c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph${{subgraph_id}}/{params.prefix}.{wildcards.chr}.subgraph_${{subgraph}}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.variant.path >> {output.merged_variant_path}
        done
        """

#vg gbwt --tags -Z t2t.grch38.58hifi.1064zmw.chrY.assembly.gbz | more
#TODO:assembly.gfa or assembly.rmtips.gfa?
rule form_chr_gbz:
    input:
        gfa = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gfa"
    output:
        chr_gbz = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gbz"
    shell:
        """
        vg gbwt --set-tag "reference_samples=CHM13" -G {input.gfa} --gbz-format -g {output.chr_gbz} --max-node 0 -p
        
        """
    
    
    
rule form_merged_gbz:
    input:
        chr_gfas = expand("c7_graph_construction/graph_merge/{prefix}.{chr}.assembly.gfa", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"], prefix = config['prefix']),
        chr_gbzs = expand("c7_graph_construction/graph_merge/{prefix}.{chr}.assembly.gbz", chr = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"], prefix = config['prefix'])
    output:
        merge_xg = f"c7_graph_construction/graph_merge/{config['prefix']}.nopath.xg",
        merge_gbwt = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbwt",
        merge_gbz = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbz"
    params:
        prefix = lambda wildcards: config['prefix']
    shell:
        """
        > c7_graph_construction/graph_merge/{params.prefix}.nopath.gfa
        for chr in chr{{{{1..22}},X,Y,M}}
        do
            grep "^S\\|^L" c7_graph_construction/graph_merge/{params.prefix}.${{chr}}.assembly.gfa >> c7_graph_construction/graph_merge/{params.prefix}.nopath.gfa
        done
        vg convert -g c7_graph_construction/graph_merge/{params.prefix}.nopath.gfa -x > {output.merge_xg}
        
        for chr in chr{{{{1..22}},X,Y,M}}
        do
            vg gbwt -Z c7_graph_construction/graph_merge/{params.prefix}.${{chr}}.assembly.gbz -o c7_graph_construction/graph_merge/{params.prefix}.${{chr}}.assembly.gbwt
        done
        vg gbwt -f c7_graph_construction/graph_merge/{params.prefix}.chr{{{{1..22}},X,Y,M}}.assembly.gbwt -o {output.merge_gbwt}
        
        vg gbwt -x {output.merge_xg} {output.merge_gbwt} --gbz-format -g {output.merge_gbz}
        """
        
# rule haplotype_path:
#     input:
#         chr_gbz = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.gbz"
#     output:
#         chr_dist = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.dist",
#         chr_ri = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.ri",
#         chr_hapl = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.hapl"
#     shell:
#         """
#         vg index -j {output.chr_dist} -t 5 {input.chr_gbz} -p
#         vg gbwt -r {output.chr_ri} -Z {input.chr_gbz}
#         vg haplotypes -t 5 --subchain-length 50000 -H {output.hapl} {input.chr_gbz} -v 1

#         """

rule haplotype_path:
    input:
        gbz = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbz"
    output:
        dist = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.dist",
        ri = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.ri",
        hapl = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.hapl"
    threads: 25
    shell:
        """
        vg index -j {output.dist} -t {threads} {input.gbz} -p
        vg gbwt -r {output.ri} -Z {input.gbz}
        vg haplotypes -t {threads} --subchain-length 50000 -H {output.hapl} {input.gbz} -v 1

        """
        
