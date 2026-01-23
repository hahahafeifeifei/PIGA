rule all_simplify_pangenome:
    input:
        lambda wildcards: expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.gfa", id = get_already_subgraph_ids(wildcards))

rule graph_bed_filtering:
    input:
        fa = get_subgraph_fa,
        gfa = get_subgraph_gfa
    output:
        high_coverage_bed = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.high_coverage.bed",
        minigraph_unaligned_bed = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.minigraph_unaligned.bed",
        merged_bed = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.merge_filter.bed",
        bed_clip_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.bed_clip.gfa",
        filter_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.gfa"
    resources:
        mem_mb=lambda wildcards, attempt: 200 * 1024 * attempt,
        runtime_hrs=lambda wildcards, attempt: 20 * attempt
    threads: 8
    shell:
        """
        sample_number=$(grep ">" {input.fa} | awk '{{split($1,a,".");if(length(a)==2) print a[1];else print a[1]"."a[2] }}' | sort -u | wc -l)
        python3 scripts/simplify_pangenome/gfa_high-coverage_bed.py {input.gfa} $[sample_number*20] 200 {output.high_coverage_bed}
        python3 scripts/simplify_pangenome/gfa_minigraph_unaligned.py {input.gfa} {output.minigraph_unaligned_bed}
        cat {output.high_coverage_bed} {output.minigraph_unaligned_bed} | bedtools sort | bedtools merge -d 1000 > {output.merged_bed}
        python3 scripts/simplify_pangenome/gfa_bed_filter.py {input.gfa} {output.merged_bed} CHM13,GRCh38,_MINIGRAPH_ {output.bed_clip_gfa}
        python3 scripts/simplify_pangenome/gfa_remove_ac0.py {output.bed_clip_gfa} {output.filter_gfa}
        """

train_sample_map = {}
with open(config["train_sample_list"]) as f:
    for line in f:
        if not line.strip():
            continue
        train_sample = line.strip().split()[0]
        truth_sample = line.strip().split()[1]
        train_sample_map[train_sample] = truth_sample

rule prepare_train_gfa:
    input:
        filter_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.gfa"
    output:
        train_raw_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.train_raw.gfa",
        train_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.train.gfa"
    params:
        all_sample = "\\t\\|".join(list(train_sample_map.keys()) + list(train_sample_map.values()))
    resources:
        mem_mb=30 * 1024
    shell:
        """
        grep -v $'{params.all_sample}\\t' {input.filter_gfa} > {output.train_raw_gfa}        
        python3 scripts/simplify_pangenome/gfa_remove_ac0.py {output.train_raw_gfa} {output.train_gfa} 
        """
        
rule process_train_node_edge:
    input:
        filter_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.gfa"
    output:
        sample_train_node = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.{{train_sample}}.train.node",
        sample_train_edge = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.{{train_sample}}.train.edge"
    params:
        truth_sample = lambda wildcards: train_sample_map.get(wildcards.train_sample)
    wildcard_constraints:
        train_sample = "|".join(list(train_sample_map.keys()))
    resources:
        mem_mb=30 * 1024
    shell:
        """
        python3 scripts/simplify_pangenome/ml_training_selection.py {input.filter_gfa} {wildcards.train_sample} {params.truth_sample} {output.sample_train_node} {output.sample_train_edge}
        [ -s {output.sample_train_node} ] || echo -e "empty\\t0" > {output.sample_train_node}
        [ -s {output.sample_train_edge} ] || echo -e "empty\\t0" > {output.sample_train_edge}
        """
        
rule merge_samples:
    input:
        sample_train_nodes = expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.{{sample}}.train.node", sample = train_sample_map.keys(), allow_missing=True),
        sample_train_edges = expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.{{sample}}.train.edge", sample = train_sample_map.keys(), allow_missing=True)
    output:
        merged_train_node = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.train.node",
        merged_train_edge = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.train.edge"
    threads: 4
    resources:
        mem_mb=60 * 1024
    shell:
        """
        sample_size=$(echo {input.sample_train_nodes} | awk '{{print NF}}')
        if [ $sample_size -gt 1 ]
        then
            csvtk -H -t join --na NONE --outer-join -f 1 {input.sample_train_nodes} | \
            awk -v OFS="\\t"  'BEGIN{{print "empty","0"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' > {output.merged_train_node}

            csvtk -H -t join --na NONE --outer-join -f 1 {input.sample_train_edges} | \
            awk -v OFS="\\t"  'BEGIN{{print "empty","0"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' > {output.merged_train_edge}
        else
            cat {input.sample_train_nodes} | \
            awk -v OFS="\\t"  'BEGIN{{print "empty","0"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' > {output.merged_train_node}
            
            cat {input.sample_train_edges} | \
            awk -v OFS="\\t"  'BEGIN{{print "empty","0"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' > {output.merged_train_edge}
        fi
        """
        
rule gnn_feature_extract:
    input:
        train_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.train.gfa",
        merged_train_node = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.train.node",
        merged_train_edge = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.train.edge"
    output:
        node_feature = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.feature",
        edge_feature = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.feature",
        node_feature_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.feature_label",
        edge_feature_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.feature_label"
    params:
        tmp_dir = "c7_graph_construction/subgraph/subgraph_{id}"
    resources:
        mem_mb=lambda wildcards, attempt: 80 * 1024 * attempt
    threads: 8
    shell:
        """
        python3 scripts/simplify_pangenome/ml_feature_extraction.py {input.train_gfa} GRCh38,CHM13 _MINIGRAPH_ {params.tmp_dir} {output.node_feature} {output.edge_feature}
        csvtk -H -t join --left-join -f 1 --na NONE {output.node_feature} {input.merged_train_node} | awk '{{if($1!="empty") print $0}}' > {output.node_feature_label}
        csvtk -H -t join --left-join -f 1 --na NONE {output.edge_feature} {input.merged_train_edge} | awk '{{if($1!="empty") print $0}}' > {output.edge_feature_label}
        """

rule statistic_prepare:
    input:
        node_feature_labels = lambda wildcards: expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.feature_label", id = get_already_subgraph_ids(wildcards)),
        edge_feature_labels = lambda wildcards: expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.feature_label", id = get_already_subgraph_ids(wildcards))
    output:
        train_node_info = f"c7_graph_construction/ml_model/{config['prefix']}.train_node.info",
        train_node_statistics = f"c7_graph_construction/ml_model/{config['prefix']}.train_node.statistics",
        train_node_dataset = f"c7_graph_construction/ml_model/{config['prefix']}.train_node.dataset",
        train_edge_info = f"c7_graph_construction/ml_model/{config['prefix']}.train_edge.info",
        train_edge_statistics = f"c7_graph_construction/ml_model/{config['prefix']}.train_edge.statistics",
        train_edge_dataset = f"c7_graph_construction/ml_model/{config['prefix']}.train_edge.dataset"
    resources:
        mem_mb=lambda wildcards, attempt: 150 * 1024 * attempt,
    threads: 8
    shell:
        """
        cat {input.node_feature_labels} | awk '{{if($NF!="NONE") print$0}}' > {output.train_node_info}
        python3 scripts/simplify_pangenome/node_prepare.py {output.train_node_info} {output.train_node_statistics} {output.train_node_dataset}

        cat {input.edge_feature_labels} | awk '{{if($NF!="NONE") print$0}}' > {output.train_edge_info}
        python3 scripts/simplify_pangenome/edge_prepare.py {output.train_edge_info} {output.train_edge_statistics} {output.train_edge_dataset}
        """

rule model_train:
    input:
        train_node_dataset = f"c7_graph_construction/ml_model/{config['prefix']}.train_node.dataset",
        train_edge_dataset = f"c7_graph_construction/ml_model/{config['prefix']}.train_edge.dataset"
    output:
        train_node_model = f"c7_graph_construction/ml_model/{config['prefix']}.node_model.pth",
        train_edge_model = f"c7_graph_construction/ml_model/{config['prefix']}.edge_model.pth",
        train_node_threshold = f"c7_graph_construction/ml_model/{config['prefix']}.node_model.threshold",
        train_edge_threshold = f"c7_graph_construction/ml_model/{config['prefix']}.edge_model.threshold"
    resources:
        mem_mb = 150*1024
    threads: 8
    shell:
        """
            python3 scripts/simplify_pangenome/model_train.py {input.train_node_dataset} {output.train_node_model} {output.train_node_threshold}
            python3 scripts/simplify_pangenome/model_train.py {input.train_edge_dataset} {output.train_edge_model} {output.train_edge_threshold}
        """

rule node_edge_inference:
    input:
        node_feature_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.feature_label",
        edge_feature_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.feature_label",
        train_node_model = f"c7_graph_construction/ml_model/{config['prefix']}.node_model.pth",
        train_edge_model = f"c7_graph_construction/ml_model/{config['prefix']}.edge_model.pth",
        train_node_threshold = f"c7_graph_construction/ml_model/{config['prefix']}.node_model.threshold",
        train_edge_threshold = f"c7_graph_construction/ml_model/{config['prefix']}.edge_model.threshold",
        train_node_statistics = f"c7_graph_construction/ml_model/{config['prefix']}.train_node.statistics",
        train_edge_statistics = f"c7_graph_construction/ml_model/{config['prefix']}.train_edge.statistics"
    output:
        model_node_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.label",
        model_edge_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.label"
    resources:
        mem_mb = 30*1024
    threads: 2
    shell:
        """
            if [ -s {input.node_feature_label} ]
            then
                node_threshold=$(cat {input.train_node_threshold})
                python3 scripts/simplify_pangenome/node_inference.py {input.node_feature_label} {input.train_node_statistics} {input.train_node_model} $node_threshold > {output.model_node_label}
            else
                > {output.model_node_label}
            fi

            if [ -s {input.edge_feature_label} ]
            then
                edge_threshold=$(cat {input.train_edge_threshold})
                python3 scripts/simplify_pangenome/edge_inference.py {input.edge_feature_label} {input.train_edge_statistics} {input.train_edge_model} $edge_threshold > {output.model_edge_label}
            else
                > {output.model_edge_label}
            fi
        """

rule gfa_ml_filter:
    input:
        node_feature_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.feature_label",
        edge_feature_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.feature_label",
        model_node_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.label",
        model_edge_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.label",
        filter_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.gfa"
    output:
        model_node_fp_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.node.fp.label",
        model_edge_fp_label = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}.subgraph_{{id}}.seqwish.smoothxg.gfaffix.filter.edge.fp.label",
        ml_filter_raw_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.raw.gfa",
        ml_filter_vg = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.vg"
    resources:
        mem_mb=lambda wildcards, attempt: 150 * 1024 * attempt
    threads: 4
    params:
        train_sample = "\\t\\|".join(list(train_sample_map.keys()))
    shell:
        """
        paste {input.node_feature_label} {input.model_node_label} | awk '{{if ($(NF-2)==0 || ($(NF-2)!=0 && $NF=="FP")) print $1}}' > {output.model_node_fp_label}
        paste {input.edge_feature_label} {input.model_edge_label} | awk '{{if ($(NF-2)==0 || ($(NF-2)!=0 && $NF=="FP")) print $1}}' > {output.model_edge_fp_label}

        python3 scripts/simplify_pangenome/gfa_node_edge_filtering.py {input.filter_gfa} {output.model_node_fp_label} {output.model_edge_fp_label} CHM13,GRCh38,_MINIGRAPH_ {output.ml_filter_raw_gfa}
        awk '{{if($1!="W" || ($2=="GRCh38" || $2=="CHM13") || $6-$5>50) print$0}}' {output.ml_filter_raw_gfa} | \
        grep -v $'{params.train_sample}\\t\\|_MINIGRAPH_\\t' | sed "s/CHM13.chr/chr/g" | \
        vg convert -g - | vg clip -d 1 -P CHM13 - | vg mod -u - > {output.ml_filter_vg}
        """

rule variant_projection:
    input:
        ml_filter_vg = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.vg",
        variant_vcf = get_merge_phased_vcf_input
    output:
        variant_project_vg = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.vg",
        variant_project_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfa",
        ref_component_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.ref_component.gfa",
        gfaffix_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.ml_filter.variant_project.gfaffix.gfa"
    resources:
        mem_mb=lambda wildcards, attempt: 150 * 1024 * attempt
    shell:
        """
        if [ $(vg paths -x {input.ml_filter_vg} -L | grep CHM13 | wc -l) -eq 0 ]
        then \
            cp {input.ml_filter_vg} {output.variant_project_vg} 
        else \
            python3 scripts/simplify_pangenome/graph_variant_projection.py \
            -g {input.ml_filter_vg} \
            -v {input.variant_vcf} \
            -r CHM13 \
            -o {output.variant_project_vg} \
            -s .snv \
            -p {wildcards.id}_snv_
        fi
        
        vg convert -f {output.variant_project_vg} > {output.variant_project_gfa}
        python3 scripts/simplify_pangenome/gfa_ref_component.py {output.variant_project_gfa} {output.ref_component_gfa} CHM13
        gfaffix {output.ref_component_gfa} -x CHM13 | vg convert -g - | vg mod -X 1024 - | vg view - > {output.gfaffix_gfa}
        """


