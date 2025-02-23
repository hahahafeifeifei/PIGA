
chr_list= [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# rule all:
#     input:
#         "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.mergeassembly.hapl"

#for a given chromosome, generate all subgraph_id.
def get_all_subgraph_assembly_gfa_files(wildcards, prefix):
    chr_subgraph_combination = checkpoints.check_chr_subgraph_combination.get().output[0]
    chr_pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, subgraph_id = line.strip().split("\t")
            if chr == wildcards.chr:
                chr_pairs.append(f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.assembly.gfa")
    return chr_pairs
        
###TODO:how to get the ordered subgraph list?
rule subgraph_gfa_merge:
    input:
        partial(get_all_subgraph_assembly_gfa_files, prefix = config['prefix']),
        subgraph_order_list = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/subgraph_order.list"
    output:
        merged_gfa = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.gfa"
    threads: 8
    shell:
        """
        awk -v chr={wildcards.chr} '{if($1==chr) print$2}' {input.subgraph_order_list} | while read subgraph_id
        do
        echo c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph${{subgraph_id}}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_${{subgraph}}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.assembly.gfa
        done > c7_graph_construction/graph_merge/{wildcards.chr}.gfa.list
        
        python3 scripts/graph-simplification/gfa_merge.py \
        c7_graph_construction/graph_merge/{wildcards.chr}.gfa.list \
        {output.merged_gfa} \
        c7_graph_construction/graph_merge/tmp.merge.{wildcards.chr} \
        {threads}
        """
#for a given chromosome, generate all subgraph_id.
def get_all_subgraph_variant_path_files(wildcards, prefix):
    chr_subgraph_combination = checkpoints.check_chr_subgraph_combination.get().output[0]
    chr_pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, subgraph_id = line.strip().split("\t")
            if chr == wildcards.chr:
                chr_pairs.append(f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.variant.path")
    return chr_pairs        
        
rule subgraph_variant_path_merge:
    input:
        partial(get_all_subgraph_variant_path_files, prefix = config['prefix']),
        subgraph_order_list = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/subgraph_order.list"
    output:
        merged_variant_path = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.variant.path"
    shell:
        """
        > {output.merged_variant_path}
        awk -v chr={wildcards.chr} '{if($1==chr) print$2}' {input.subgraph_order_list} | while read subgraph_id
        do
        cat c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph${{subgraph_id}}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_${{subgraph}}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.variant.path >> {output.merged_variant_path}
        done
        """

#vg gbwt --tags -Z t2t.grch38.58hifi.1064zmw.chrY.assembly.gbz | more
#TODO:assembly.gfa or assembly.rmtips.gfa?
rule form_chr_gbz:
    input:
        gfa = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.gfa"
    output:
        chr_gbz = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.gbz"
    shell:
        """
        vg gbwt --set-tag "reference_samples=CHM13" -G {input.gfa} --gbz-format -g {output.chr_gbz} --max-node 0 -p
        
        """
    
    
    
#TODO: where do the gbz files come from?
rule form_merged_gbz:
    input:
        chr_gfas = expand("c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.gfa", chr=chr_list),
        chr_gbzs = expand("c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.gbz", chr=chr_list)
    output:
        merge_xg = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.nopath.xg",
        merge_gbwt = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.gbwt",
        merge_gbz = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.gbz"
    shell:
        """
        > c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.nopath.gfa
        for chr in chr{{{{1..22}},X,Y,M}}
        do
            grep "^S\\|^L" c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.${{chr}}.assembly.gfa >> c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.nopath.gfa
        done
        vg convert -g c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.nopath.gfa -x > {output.merge_xg}
        
        for chr in chr{{{{1..22}},X,Y,M}}
        do
            vg gbwt -Z c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.${{chr}}.assembly.gbz -o c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.${{chr}}.assembly.gbwt
        done
        vg gbwt -f c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.chr{{{{1..22}},X,Y,M}}.assembly.gbwt -o {output.merge_gbwt}
        
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
        gbz = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.gbz"
    output:
        dist = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.dist",
        ri = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.ri",
        hapl = "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.hapl"
    threads: 25
    shell:
        """
        vg index -j {output.dist} -t 25 {input.gbz} -p
        vg gbwt -r {output.ri} -Z {input.gbz}
        vg haplotypes -t 25 --subchain-length 50000 -H {output.hapl} {input.gbz} -v 1

        """
        
