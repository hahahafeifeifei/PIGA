
chr_list= [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

id_list = get_already_subgraph_ids()
wildcard_constraints:
    id = "|".join(id_list)

rule all_merge_pangenome:
    input:
        f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.hapl",
        f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbz",
        expand(f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gfa", chr = chr_list),
        expand(f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.variant.path", chr = chr_list)


rule subgraph_feature:
    input:
        gfaffix_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.gfa"
    output:
        node_count = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.node.count",
        chr_subgraph = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.chr_subgraph.list"
    resources:
        mem_mb = 60*1024
    shell:
        """
        grep ^S {output.gfaffix_gfa} | wc -l | awk -v id={wildcards.id} '{{print id"\\t"$1}}' > {output.node_count}
        chr=$(grep CHM13 {output.gfaffix_gfa} | head -1 | awk '{{split($4,a,".");print a[length(a)]}}')
        echo -e $chr"\\t"{wildcards.id} > {output.chr_subgraph}
        """

rule feature_merge:
    input:
        counts = expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.node.count", id=id_list),
        chr_subgraphs = expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.chr_subgraph.list", id=id_list)
    output:
        node_counts = f"c7_graph_construction/{config['prefix']}.node.count",
        node_sum_counts = f"c7_graph_construction/{config['prefix']}.node.sum.count",
        chr_subgraph_list = f"c7_graph_construction/{config['prefix']}.chr_subgraph.list"
    resources:
        mem_mb = 60*1024
    shell:
        """
        cat {input.counts} > {output.node_counts}
        awk 'BEGIN{{sum=0}}{{print $1"\\t"sum;sum+=$2}}' {output.node_counts} > {output.node_sum_counts}
        cat {input.chr_subgraphs} > {output.chr_subgraph_list}
        """

rule node_ids:
    input:
        gfaffix_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.gfa",
        node_sum_counts = f"c7_graph_construction/{config['prefix']}.node.sum.count"
    output:
        ids_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.ids.gfa",
        ids_assembly_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.ids.assembly.gfa",
        ids_variant_path = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.ids.variant.path"
    resources:
        mem_mb = 60*1024,
        runtime_hrs=10
    shell:
        """
        index=$(awk -v id={wildcards.id} '{{if($1==id) print $2}}' {input.node_sum_counts})
        python3 scripts/simplify_pangenome/gfa_ids.py {input.gfaffix_gfa} {output.ids_gfa} $index
        grep -v snv {output.ids_gfa} > {output.ids_assembly_gfa}
        grep snv {output.ids_gfa} > {output.ids_variant_path}
        """


rule chr_merge:
    input:
        chr_subgraph_list = f"c7_graph_construction/{config['prefix']}.chr_subgraph.list"
    output:
        gfa_list = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.gfa.list"
        merged_gfa = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gfa",
        merged_variant_path = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.variant.path"
    threads: 8
    wildcard_constraints:
        chr = "|".join(chr_list)
    params:
        prefix = config['prefix']
    shell:
        """
        awk -v chr={wildcards.chr} '{{if($1==chr) print$2}}' {input.chr_subgraph_list} | while read id
        do
            echo c7_graph_construction/subgraph/subgraph_${id}/{params.prefix}_subgraph_${id}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.ids.assembly.gfa
        done > {input.gfa_list}
        
        python3 scripts/simplify_pangenome/gfa_merge.py {input.gfa_list} {output.merged_gfa} \
        c7_graph_construction/graph_merge/tmp.merge.{wildcards.chr} CHM13 {threads}

        awk -v chr={wildcards.chr} '{{if($1==chr) print$2}}' {input.chr_subgraph_list} | while read id
        do
            cat c7_graph_construction/subgraph/subgraph_${id}/{params.prefix}_subgraph_${id}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.ids.variant.path
        done > {output.merged_variant_path}
        """


rule chr_gbz:
    input:
        chr_gfa = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gfa"
    output:
        chr_gbz = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gbz",
        chr_gbwt = f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gbwt"
    threads: 4
    shell:
        """
        vg gbwt --set-tag "reference_samples=CHM13" -G {input.chr_gfa} --gbz-format -g {output.chr_gbz} --max-node 0
        vg gbwt -Z {output.chr_gbz} -o {output.chr_gbwt}
        """
    
rule merged_gbz:
    input:
        chr_gfas = expand(f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gfa", chr = chr_list),
        chr_gbzs = expand(f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gbz", chr = chr_list),
        chr_gbwts = expand(f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.assembly.gbwt", chr = chr_list)
    output:
        merge_gfa = f"c7_graph_construction/graph_merge/{config['prefix']}.nopath.gfa",
        merge_xg = f"c7_graph_construction/graph_merge/{config['prefix']}.nopath.xg",
        merge_gbwt = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbwt",
        merge_gbz = f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbz"
    threads: 25
    shell:
        """
        grep "^S\\|^L" {input.chr_gfas} > {output.merge_gfa}
        vg convert -g {output.merge_gfa} -x > {output.merge_xg}
        vg gbwt -f {input.chr_gbwts} -o {output.merge_gbwt}
        vg gbwt -x {output.merge_xg} {output.merge_gbwt} --gbz-format -g {output.merge_gbz}
        """

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
        vg index -j {output.dist} -t {threads} {input.gbz} --no-nested-distance
        vg gbwt -r {output.ri} -Z {input.gbz}
        vg haplotypes -t {threads} --subchain-length 50000 --linear-structure -H {output.hapl} {input.gbz} -v 1
        """
        
