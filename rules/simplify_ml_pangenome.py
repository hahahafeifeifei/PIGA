



            
def get_all_subgraph_assembly_gfa_files(wildcards, prefix):
    chr_subgraph_combination = checkpoints.check_chr_subgraph_combination.get().output[0]
    pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, subgraph_id = line.strip().split("\t")
            pairs.append(f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.assembly.gfa")
    return pairs

def get_all_subgraph_variant_path_files(wildcards, prefix):
    chr_subgraph_combination = checkpoints.check_chr_subgraph_combination.get().output[0]
    pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, subgraph_id = line.strip().split("\t")
            pairs.append(f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.variant.path")
    return pairs 

# rule all:
#         input:
#             partial(get_all_subgraph_assembly_gfa_files, config['prefix']),
#             partial(get_all_subgraph_variant_path_files, config['prefix'])
            
            
rule prepare_training_set:
    input:
        linear_gfaffix_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfa"
    output:
        training_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.training.gfa"
    params:
        exclude_samples = "|".join(config['training_samples'] + config['test_samples'])
    shell:
        """
        grep -vE '{params.exclude_samples}\\t' \
        {input.linear_gfaffix_gfa} \
        > c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{wildcards.subgraph_id}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.remove.gfa
        
        python3 scripts/graph-simplification/gfa_remove_ac0.py c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.remove.gfa {output.training_gfa}
        """
        
rule process_validation_sample:
    input:
        linear_gfaffix_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfa"
    output:
        sample_training_node = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.training.node",
        sample_training_edge = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.training.edge"
    shell:
        """
        python3 scripts/graph-simplification/gnn_training_selection.py {input.gfa} {wildcards.sample}_low {wildcards.sample} {output.sample_training_node} {output.training_edge}
        [ -s {output.sample_training_node} ] || echo -e "empty\\t0" > {output.sample_training_node}
        [ -s {output.sample_training_edge} ] || echo -e "empty\\t0" > {output.sample_training_edge}
        """
        
rule merge_samples:
    input:
        sample_training_nodes = expand("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.training.node",sample=config['validate_samples'], prefix=config['prefix'], allow_missing=True),
        sample_training_edges = expand("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.training.edge",sample=config['validate_samples'], prefix=config['prefix'], allow_missing=True)
    output:
        merged_training_node = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.merge.training.node",
        merged_training_edge = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.merge.training.edge",
        merged_all_node = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.merge.all.node",
        merged_all_edge = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.merge.all.edge"
        
    shell:
        """
        csvtk -H -t join --na NONE --outer-join -f 1 {input.sample_training_nodes} | \
        awk 'BEGIN{{OFS="\\t"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} \
            if((tp>0 && fp==0) || (tp==0 && fp>0)) print $0 }}' | cut -f 1,2,3,4,5,7 | grep -v empty | \
        awk 'BEGIN{{OFS="\\t"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} \
            if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' > {output.merged_training_node}
        
        csvtk -H -t join --na NONE --outer-join -f 1 {input.sample_training_edges} | \
        awk 'BEGIN{{OFS="\\t"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} \
            if((tp>0 && fp==0) || (tp==0 && fp>0)) print $0 }}' | cut -f 1,2,3,4,5,7 | grep -v empty | \
        awk 'BEGIN{{OFS="\\t"}} {{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} \
            if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' > {output.merged_training_edge}
            
        csvtk -H -t join --na NONE --outer-join {input.sample_training_nodes} | \
        grep -v empty | \
        awk -v OFS='\\t' '{{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' \
        > {output.merged_all_node}
        
        csvtk -H -t join --na NONE --outer-join {input.sample_training_edges} | \
        grep -v empty | \
        awk -v OFS='\\t' '{{tp=0;fp=0;for(i=2;i<=NF;i++) {{if($i=="TP")tp+=1;if($i=="FP")fp+=1}} if(tp>0 && fp==0) print $1,"1"; if(tp==0 && fp>0) print $1,"0" }}' \
        > {output.merged_all_edge}
        """
        
rule generate_validation_files:
    input:
        merged_training_node = rules.merge_samples.output.merged_training_node,
        merged_training_edge = rules.merge_samples.output.merged_training_edge,
        merged_all_node = rules.merge_samples.output.merged_all_node,
        merged_all_edge = rules.merge_samples.output.merged_all_edge
    output:
        validation_node = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.merge.validation.node"),
        validation_edge = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.merge.validation.edge")
    shell:
        """
        csvtk -H -t join --na NONE --left-join \
        {input.merged_all_node} \
        {input.merged_training_node} | \
        awk -v OFS='\\t' '{{if(($2==0 && $3=="NONE") || ($2==1 && $3=="NONE"))print$1,$2}}' > {output.validation_node}

        csvtk -H -t join --na NONE --left-join \
        {input.merged_all_edge} \
        {input.merged_training_edge} | \
        awk -v OFS='\\t' '{{if(($2==0 && $3=="NONE") || ($2==1 && $3=="NONE"))print$1,$2}}' > {output.validation_edge}
        
        [ -s {input.merged_training_node} ] || echo -e "empty\\t0" > {input.merged_training_node}
        [ -s {input.merged_training_edge} ] || echo -e "empty\\t0" > {input.merged_training_edge}
        [ -s {output.validation_node} ] || echo -e "empty\\t0" > {output.validation_node}
        [ -s {output.validation_edge} ] || echo -e "empty\\t0" > {output.validation_edge}
        
        """

rule generate_special_validation:
    input:
        sample_training_node = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.training.node",
        sample_training_edge = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.training.edge"
    output:
        sample_validation_node = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.validation.node",
        sample_validation_edge = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.validation.edge"
    shell:
        """
        awk -v OFS='\\t' '{{if($2=="TP")print $1,1;else print $1,0}}' {input.sample_training_node} > {output.sample_validation_node}
        awk -v OFS='\\t' '{{if($2=="TP")print $1,1;else print $1,0}}' {input.sample_training_edge} > {output.sample_validation_edge}
        """
        
rule gnn_feature_extract:
    input:
        training_gfa = rules.prepare_training_set.output.training_gfa,
    output:
        node_features = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.training.node",
        edge_features = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.training.edge"
    shell:
        """
        python3 {script_dir}/gnn_feature_extraction.py \
        {input.training_gfa} \
        {output.node_features} \
        {output.edge_features}
        
        """

#TODO: currently I assume only two samples in special_validate_samples. what if more than 2?
rule finalize_labels:
    input:
        node_features = rules.gnn_feature_extract.output.node_features,
        merged_training_node = rules.merge_samples.output.merged_training_node,
        merged_validation_node = rules.generate_validation_files.output.validation_node,
        test_validation_node = expand("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.validation.node", sample = config['test_samples'], prefix=config['prefix'], allow_missing=True),
        edge_features = rules.gnn_feature_extract.output.edge_features,
        merged_training_edge = rules.merge_samples.output.merged_training_edge,
        merged_validation_edge = rules.generate_validation_files.output.validation_edge,
        test_validation_edge = expand("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.{sample}.validation.edge", sample = config['test_samples'], prefix=config['prefix'], allow_missing=True)
    output:
        node_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.training.node.label",
        edge_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.training.edge.label"
    shell:
        """
        csvtk -H -t join --left-join -f 1 --na NONE \
        {input.node_features} \
        {input.merged_training_node} | \
        csvtk -H -t join --left-join -f 1 --na NONE - {input.test_validation_node}[0] \
        csvtk -H -t join --left-join -f 1 --na NONE - {input.test_validation_node}[1] \
        > {output.node_label}
        
        csvtk -H -t join --left-join -f 1 --na NONE \
        {input.edge_features} \
        {input.merged_training_edge} | \
        csvtk -H -t join --left-join -f 1 --na NONE - {input.test_validation_edge}[0] \
        csvtk -H -t join --left-join -f 1 --na NONE - {input.test_validation_edge}[1] \
        > {output.edge_label}
        """
        
### training part

#TODO: column 83,84 should be adjusted.
rule statistic_prepare:
    input:
        node_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.training.node.label",
        edge_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.training.edge.label"
    output:
        train_node_info = "c7_graph_construction/chr_mc/{chr}/subgraph/train_node.info",
        train_node_statistics = "c7_graph_construction/chr_mc/{chr}/subgraph/train_node.statistics",
        train_node_dataset = "c7_graph_construction/chr_mc/{chr}/subgraph/train_node.dataset",
        train_edge_info = "c7_graph_construction/chr_mc/{chr}/subgraph/train_edge.info",
        train_edge_statistics = "c7_graph_construction/chr_mc/{chr}/subgraph/train_edge.statistics",
        train_edge_dataset = "c7_graph_construction/chr_mc/{chr}/subgraph/train_edge.dataset"
    shell:
        """
        cat chr_mc/{wildcards.chr}/subgraph/*/*linearize.training.node.label | awk '{{if($83!="NONE" || $84!="NONE") print$0}}'> {output.train_node_info}
        python3 scripts/graph-simplification/node_prepare.py \
        {output.train_node_info} \
        {output.train_node_statistics} \
        {output.train_node_dataset}
        
        cat chr_mc/{wildcards.chr}/subgraph/*/*linearize.training.edge.label | awk '{{if($83!="NONE" || $84!="NONE") print$0}}'> {output.train_edge_info}
        python3 scripts/graph-simplification/edge_prepare.py \
        {output.train_edge_info} \
        {output.train_edge_statistics} \
        {output.train_edge_dataset}
        """

#TODO: model_dir path? How to get the model?
rule node_edge_inference:
    input:
        node_label = rules.finalize_labels.output.node_label,
        edge_label = rules.finalize_labels.output.edge_label,
        node_statistics = rules.statistic_prepare.output.train_node_statistics,
        edge_statistics = rules.statistic_prepare.output.train_edge_statistics
        
    output:
        TVR90_node_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.inference.TVR90.node.label",
        TVR90_edge_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.inference.TVR90.edge.label"
    params:
        chr = "{chr}",
        model_dir = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/model-train/version2"
    run:
        if params.chr == "chrX":
            node_inference_script = f"{params.model_dir}/node_X/node_inference.py",
            edge_inference_script = f"{params.model_dir}/edge_X/edge_inference.py",
            node_model = f"{params.model_dir}/node_X/node.X.pth",
            edge_model = f"{params.model_dir}/edge_X/edge.X.pth",
            node_threshold = 0.0371781,
            edge_threshold = 0.04142026
        elif params.chr == "chrY":
            node_inference_script = f"{params.model_dir}/node_Y/node_inference.py",
            edge_inference_script = f"{params.model_dir}/edge_Y/edge_inference.py",
            node_model = f"{params.model_dir}/node_Y/node.Y.pth",
            edge_model = f"{params.model_dir}/edge_Y/edge.Y.pth",
            node_threshold = 0.08506753,
            edge_threshold = 0.06190881
        else:
            node_inference_script = f"{params.model_dir}/node_auto/node_inference.py",
            edge_inference_script = f"{params.model_dir}/edge_auto/edge_inference.py",
            node_model = f"{params.model_dir}/node_auto/node.auto.pth",
            edge_model = f"{params.model_dir}/edge_auto/edge.auto.pth",
            node_threshold = 0.04559464,
            edge_threshold = 0.04971138
        shell("python3 {node_inference_script} {input.node_label} {input.node_statistics} {node.model} {node_threshold} > {output.TVR90_node_label}")
        shell("python3 {edge_inference_script} {input.edge_label} {input.edge_statistics} {edge.model} {edge_threshold} > {output.TVR90_edge_label}")


rule FP_node_edge_prepare:
    input:
        node_label = rules.finalize_labels.output.node_label,
        TVR90_node_label = rules.node_edge_inference.output.TVR90_node_label,
        edge_label = rules.finalize_labels.output.edge_label,
        TVR90_edge_label = rules.node_edge_inference.output.TVR90_edge_label
    output:
        TVR90_FP_node_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.inference.TVR90.FP.node",
        TVR90_FP_edge_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.inference.TVR90.FP.edge"
    shell:
        """
        paste {input.node_label} {input.TVR90_node_label} | awk '{{if($83==0 || $84==0) print$87;else {{if($88=="FP" && $83!=1 && $84!=1) print$87}} }}' > {output.TVR90_FP_node_label}
        paste {input.edge_label} {input.TVR90_edge_label} | awk '{{if($78==0 || $79==0) print$82;else {{if($83=="FP" && $78!=1 && $79!=1) print$82}} }}' > {output.TVR90_FP_edge_label}
        """

#TODO:edit the script to removed the parameters hifi_samples_list.
rule gfa_node_edge_TVR90_filter:
    input:
        linear_gfaffix_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfa",
        TVR90_FP_node_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.inference.TVR90.FP.node",
        TVR90_FP_edge_label = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.inference.TVR90.FP.edge"
    output:
        TVR90_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.gfa",
        TVR90_filter_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.filter.gfa",
        TVR90_rmac0_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.rmac0.gfa"
    params:
        # hifi_samples_list = ",".join(config['hifi_samples']),
        validate_low_samples_list = "|".join(f"{sample}_low\t" for sample in config['validate_samples'])
    shell:
        """
        python3 scripts/graph-simplification/gfa_node_edge_filtering.py \
        {input.linear_gfaffix_gfa} \
        {input.TVR90_FP_node_label} \
        {input.TVR90_FP_edge_label} \
        {params.hifi_samples_list} \
        {output.TVR90_gfa}
        
        awk '{{if($1!="W" || ($2=="GRCh38" || $2=="CHM13") || $6-$5>50) print$0}}' \
        {output.TVR90_gfa} | grep -v '{params.validate_low_samples_list}|_MINIGRAPH_\t' > {output.TVR90_filter_gfa}
        
        python3 scripts/graph-simplification/gfa_remove_ac0.py \
        {output.TVR90_filter_gfa} \
        {output.TVR90_filter_rmac0_gfa}
        """
        


#TODO: can be written by pipe.    
# params:
#     hifi_samples_dashP_command = "|".join((f"-P {sample} " for sample in config['hifi_samples']))
rule TVR90_vg_clip:
    input:
        TVR90_rmac0_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.rmac0.gfa"
    output:
        vg = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.rmac0.vg"),
        vg_clip = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.rmtips.vg"),
        vg_unchop = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.vg"),
        gfa_unchop = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.gfa"
    shell:
        """
        vg convert -g {input.TVR90_rmac0_gfa} > {output.vg}
        vg clip {output.vg} > {output.vg_clip}
        vg mod -u {output.vg_clip} > {output.vg_unchop}
        vg view {output.vg_unchop} > {output.gfa_unchop}
        
        """
        
#can iterate multiple times.
rule TVR90_snarls_filter:
    input:
        gfa = lambda wildcards: f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_{wildcards.subgraph}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{int(wildcards.i)-1}_filter.gfa",
        stat = lambda wildcards: f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_{wildcards.subgraph}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{int(wildcards.i)-1}.stat"
    output:
        out_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{i}_filter.gfa",
        out_stat = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{i}.stat",
        filter = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{i}_filter.bed",
        txt = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{i}.txt"),
        region = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{i}_region.bed"),
        temp1_gfa = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{i}_tmp1.gfa")
    shell:
        """
        if [ $(awk '{{if($7>20000 && $8>20000) print$0}}' {input.stat} | wc -l) -eq 1 ]; then \
            cp {input.gfa} {output.out_gfa}
            cp {input.stat} {output.out_stat}
        else \
            vg snarls -a -T {input.gfa} | vg view -R - > {output.txt}
            
            python3 scripts/graph-simplification/gfa_snarl_bed.py \
                {input.gfa} {input.stat} {output.txt} \
                20000 100 {output.region}
                
            bedtools sort -i {output.region} | bedtools merge -i - -d 50 > {output.filter}
            
            python3 scripts/graph-simplification/gfa_bed_filter.py {input.gfa} {output.filter} CHM13,GRCh38 {output.temp1_gfa}
            
            python3 scripts/graph-simplification/gfa_remove_ac0.py {output.temp1_gfa} {output.gfa}
            
            vg stats -R {output.out_gfa} > {output.out_stat}
        fi
        """

rule component_extract:
    input:
        snarls20_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls20_filter.gfa"
    output:
        out_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls20_filter.component.gfa",
        out_stat = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls20_filter.component.snarls.stat"
    shell:
        """
        python3 scripts/graph-simplification/gfa_ref_component.py \
        {input.snarls20_gfa} \
        {output.out_gfa} \
        CHM13,GRCh38
        
        vg stats -R {output.out_gfa} > {output.out_stat}
        """
        
#TODO: where do the snarls_filter.ref_node.gfa come from?
rule ref_node_filter:
    input:
        bed = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls20_filter.bed",
        ref_node_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls_filter.ref_node.gfa"
        
    output:
        bed = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls_filter.bed",
        gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls_filter.gfa"
    shell:
        """
        cat c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls{1..20}_filter.bed | bedtools sort -i - | bedtools merge -i - > {output.bed}
        awk '{{if($1!="W" || ($2=="GRCh38" || $2=="CHM13") || $6-$5>50) print$0}}' {ref_node_gfa} > {output.gfa}
        
        """

#TODO: path of variants should be adjusted.
rule variant_projection:
    input:
        gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls20_filter.component.gfa",
        variants = "/storage/yangjianLab/wangyifei/project/01.CKCG/11.Consensus_variant/CHM13/03.shapeit/bench_set/merge_vcf/CKCG.consensus.whatshap.shapeit4.mutiallele.sample_match.bench_set.vcf.gz"
    output:
        vg = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls20_filter.component.vg",
        clip_vg = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.unchop.snarls20_filter.component.rmtips.vg",
        variant_project_vg = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.vg",
        variant_project_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfa"
    shell:
        """
        awk '{{if($1!="W" || ($2=="GRCh38" || $2=="CHM13") || $6-$5>50) print$0}}' {input.gfa} | vg convert -g - > {output.vg}
        
        vg clip \
        {output.vg} \
        > {output.clip_vg}
        
        if [ {wildcards.chr} == "chrY" ]
        then \
            cp {output.clip_vg} {output.variant_project_vg} 
        else \
            python3 scripts/graph-simplification/graph_variant_projection.py \
            -g {output.clip_vg} \
            -v {input.variants} \
            -r CHM13 \
            -o {output.variant_project_vg} \
            -s .snv \
            -p {wildcards.chr}_{wildcards.subgraph_id}_snv_
        fi
        
        vg convert -fW {output.variant_project_vg} > {output.variant_project_gfa}
        """
        
rule variant_project_gfa_gfaffix:
    input:
        variant_project_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfa"
    output:
        gfaffix_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.gfa",
        gfaffix_info = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.info"
    shell:
        """
        gfaffix {input.variant_project_gfa} -o {output.gfaffix_gfa} > {output.gfaffix_info}
        
        """
        
rule variant_project_gfa_chop:
    input:
        gfaffix_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.gfa"
    output:
        gfaffix_vg = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.vg"),
        unchop_vg = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.unchop.vg"),
        chop_vg = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.vg"),
        chop_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.gfa"
    shell:
        """
        vg convert -g {output.gfaffix_gfa} > {output.gfaffix_vg}
        
        vg mod -u {output.gfaffix_vg} > {output.unchop_vg}
        vg mod -X 1024 {output.unchop_vg} > {output.chop_vg}
        
        vg view {output.chop_vg} > c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}/t2t.grch38.58hifi.1064zmw.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.gfa
        """

#TODO: how to generate the node_info for each subgraph?
rule variant_project_gfa_ids:
    input:
        node_list = "/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/subgraph_order.chop.node.sum.list",
        chop_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.gfa"
    output:
        ids_gfa = temp("c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.gfa"),
        ids_assembly_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.assembly.gfa",
        ids_variant_path = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/t2t.grch38.58hifi.1064zmw.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.variant.path"
    shell:
        """
        node_index=$(awk -v chr=$chr -v subgraph=$subgraph '{{if($1==chr && $2==subgraph) print$3}}' {input.node_list})
        
        python3 scripts/graph-simplification/gfa_ids.py \
        {input.chop_gfa} \
        {output.ids_gfa} \
        $node_index
        
        grep -v snv {output.ids_gfa} > {output.ids_assembly_gfa}
        grep snv {output.ids_gfa} > {output.ids_variant_path}
        """        

# rule xxx:
#     input:
#     output:
#     shell:
#         """
        
#         """  