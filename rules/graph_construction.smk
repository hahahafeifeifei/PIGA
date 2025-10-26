chr_list = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


#prefix = '1kcp'
            
def get_linear_gfaffix_gfa_files(wildcards, prefix):
    if "subgraph_list" in config:
        chr_subgraph_combination_file = config["subgraph_list"]
    else:
        chr_subgraph_combination_file = checkpoints.prepare_subgraph_list.get().output[0]
    pairs = []
    with open(chr_subgraph_combination_file) as f:
            for line in f:
                chr_id, subgraph_id = line.strip().split("\t")
                pairs.append(f"c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}/{prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfa")
    return pairs
            
rule all_graph_construction:
    input:
        partial(get_linear_gfaffix_gfa_files, prefix = config['prefix'])
    
 

    
    
rule subgraph_seqwish:
    input:
        paf = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.paf"
        
    output:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.gfa"
    threads: 1
    params:
        prefix = lambda wildcards: config['prefix'],
        chr_subgraph_dir = "c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}"
    shell:
        """
        cat {params.chr_subgraph_dir}/fa/*.fasta > {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.fasta
        seqwish -P -t {threads} \
            -s {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.fasta \
            -p {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.paf \
            -g {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.origin.gfa
        
        awk -v OFS='\t' '{{if(substr($1,1,1)=="P") {{split($2,a,"id=");split(a[2],b,"|");if(b[1]=="CHM13" || b[1]=="GRCh38" || b[1]=="_MINIGRAPH_") name=b[1]"#"b[2];else {{split(b[1],c,".");name=c[1]"#"c[2]-1"#"b[2]"#0"}}  print $1,name,$3,$4 }}else  print$0  }}' {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.origin.gfa > {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.gfa
        
        rm {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.fasta
        rm {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.origin.gfa
        """

        
rule subgraph_smoothxg_1:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.gfa"
    output:
        smoothxg_1_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg1.gfa"
    threads: 1
    params:
        chr_subgraph_dir = "c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}"
    shell:
        """
        sample_number=$(ls -lh {params.chr_subgraph_dir}/fa/*fasta | awk '{{if($5!=0)print $0}}' | wc -l)
        mkdir {params.chr_subgraph_dir}/smoothxg_tmp
        smoothxg -t {threads} \
            -g {input.gfa} \
            -r ${{sample_number}} \
            --base {params.chr_subgraph_dir}/smoothxg_tmp \
            --chop-to 100 -I 0.98 -R 0 -j 0 -e 0 \
            -l 1400 \
            -q 2800 \
            -w $[1400*${{sample_number}}] \
            -p 1,19,39,3,81,1 -O 0.001 \
            -Y $[100*${{sample_number}}] -d 0 -D 0 -V -c 200M -W 1 \
            -o {output.smoothxg_1_gfa}
        """
        
rule subgraph_smoothxg_2:
    input:
        smoothxg_1_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg1.gfa"
        
    output:
        smoothxg_2_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfa"
    threads: 1
    params:
        chr_subgraph_dir = "c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}"
    shell:
        """
        smoothxg -t {threads} \
            -g {input.smoothxg_1_gfa} \
            -r ${{sample_number}} \
            --base smoothxg_tmp \
            --chop-to 100 -I 0.98 -R 0 -j 0 -e 0 \
            -l 1800 \
            -q 3600 \
            -w $[1800*${{sample_number}}] \
            -p 1,19,39,3,81,1 -O 0.001 \
            -Y $[100*${{sample_number}}] -d 0 -D 0 -V -c 200M -W 1 \
            -o {output.smoothxg_2_gfa}
        """
        
        
# rule subgraph_smoothxg_3:
#     input:
#         smoothxg_2_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfa"
        
#     output:
#         smoothxg_3_gfa = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg3.gfa"
#     params:
#         chr_subgraph_dir = "c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}"
#     shell:
#         """
#         smoothxg -t {threads} \
#         -g {input.smoothxg_2_gfa} \
#         -r ${{sample_number}} \
#         --base smoothxg_tmp \
#         --chop-to 100 -I 0.98 -R 0 -j 0 -e 0 \
#         -l 2200 \
#         -q 4400 \
#         -w $[2200*${{sample_number}}] \
#         -p 1,19,39,3,81,1 -O 0.001 \
#         -Y $[100*${{sample_number}}] -d 0 -D 0 -V -c 200M -W 1 \
#         -o {output.smoothxg_3_gfa}
#         """
        
        
rule subgraph_gfaffix:
    input:
        smoothxg_2_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfa" 
    output:
        gfaffix_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.gfa",
        gfaffix_info = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.info"
    shell:
        """
        gfaffix \
            {input.smoothxg_2_gfa} \
            -o {output.gfaffix_gfa} \
            > {output.gfaffix_info}
        
        """

rule gfa_format_norm:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.gfa"    
    output:
        norm_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.norm.gfa"
    params:
        prefix = lambda wildcards: config['prefix']
    shell:
        """
        vg convert -gf {input.gfa} > c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.W.gfa
        
        awk -v OFS='\t' '{{ if(substr($0,1,1)=="W" && $2!="_MINIGRAPH_") {{split($4,a,":");split(a[2],b,"-");print$1,$2,$3,a[1],$5+b[1]-1,$6+b[1]-1,$7}} else print$0 }}' \
        c7_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.W.gfa \
        > {output.norm_gfa}
        
        """
        
rule record_gfa_cycle_region:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.norm.gfa"
    output:
        cycle_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.norm.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/cycle_region.py \
            {input.gfa} \
            {output.cycle_region}
        """
        

rule record_gfa_high_cov_region:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.norm.gfa"
    output:
        high_cov_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.high_coverage.bed"
    shell:
        """
        python3 scripts/graph-simplification/gfa_high-coverage_bed.py \
            {input.gfa} \
            45200 \
            200 \
            {output.high_cov_region}
        """        

rule record_gfa_minigraph_unaligned_region:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.norm.gfa"
    output:
        minigraph_unaligned_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.minigraph_unaligned.bed"
    shell:
        """
        python3 scripts/graph-simplification/gfa_minigraph_unaligned.py \
            {input.gfa} \
            {output.minigraph_unaligned_region}
        """

rule gfa_regions_merge:
    input:
        minigraph_unaligned_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.minigraph_unaligned.bed",
        high_cov_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.high_coverage.bed"
    output:
        merged_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.merge.bed"
    shell:
        """
        cat {input.minigraph_unaligned_region} {input.high_cov_region} | bedtools sort | bedtools merge -d 1000 > {output.merged_region}
        
        """        
rule gfa_filter:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.norm.gfa",
        merged_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.merge.bed"
    output:
        filtered_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.clip.gfa" 
    shell:
        """
        python3 scripts/graph-simplification/gfa_bed_filter.py \
            {input.gfa} \
            {input.merged_region} \
            CHM13,GRCh38,_MINIGRAPH_ \
            {output.filtered_gfa}
        """        
        
        
rule record_filtered_gfa_cycle_region:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.clip.gfa"
    output:
        cycle_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.clip.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/cycle_region.py \
            {input.gfa} \
            {output.cycle_region}
        """

rule prepare_ref_file_samples_file:
    output:
        ref_file = "c7_graph_construction/ref.sample",
        sample_file = "c7_graph_construction/all.sample"
    run:
        with open(output.ref_file, 'w') as f:
            f.write('CHM13\nGRCh38\n')
        with open(output.sample_file, 'w') as f:
            f.write('\n'.join(config['samples']) + '\n')


rule gfa_linear_1:
    input:
        gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.clip.gfa",
        ref_file = "c7_graph_construction/ref.sample",
        sample_file = "c7_graph_construction/all.sample"
    output:
        linear_1_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize1.gfa",
        linear_1_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize1.region",
        linear_1_cycle_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize1.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub_muti_match.py \
            {input.gfa} 5000 50 {input.ref_file} {input.sample_file} \
            {output.linear_1_gfa} \
            {output.linear_1_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
            {output.linear_1_gfa} \
            {output.linear_1_cycle_region}
        """
        
rule gfa_linear_2:
    input:
        linear_1_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize1.gfa",
        ref_file = "c7_graph_construction/ref.sample",
        sample_file = "c7_graph_construction/all.sample"
    output:
        linear_2_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize2.gfa",
        linear_2_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize2.region",
        linear_2_cycle_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize2.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub.py \
            {input.linear_1_gfa} \
            5000 50 10000 {input.ref_file} {input.sample_file} \
            {output.linear_2_gfa} \
            {output.linear_2_cycle_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
            {output.linear_2_gfa} \
            {output.linear_2_cycle_region}
        """
        
rule gfa_linear_3:
    input:
        linear_2_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize2.gfa",
        ref_file = "c7_graph_construction/ref.sample",
        sample_file = "c7_graph_construction/all.sample"
    output:
        linear_3_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize3.gfa",
        linear_3_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize3.region",
        linear_3_cycle_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize3.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub.py \
            {input.linear_2_gfa} \
            7500 50 5000 {input.ref_file} {input.sample_file} \
            {output.linear_3_gfa} \
            {output.linear_3_cycle_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
            {output.linear_3_gfa} \
            {output.linear_3_cycle_region}
        """
        
rule gfa_linear_4:
    input:
        linear_3_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize3.gfa",
        ref_file = "c7_graph_construction/ref.sample",
        sample_file = "c7_graph_construction/all.sample"
    output:
        linear_4_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize4.gfa",
        linear_4_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize4.region",
        linear_4_cycle_region = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize4.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub.py \
            {input.linear_3_gfa} \
            9000 50 2000 {input.ref_file} {input.sample_file} \
            {output.linear_4_gfa} \
            {output.linear_4_cycle_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
            {output.linear_4_gfa} \
            {output.linear_4_cycle_region}
        """

        
        
        
rule gfa_linear_rmac0:
    input:
        linear_4_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize4.gfa"
    output:
        linear_rmac0_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize.rmac0.gfa"
    shell:
        """
        python3 scripts/graph-simplification/gfa_remove_ac0.py \
            {input.linear_4_gfa} \
            {output.linear_rmac0_gfa}
        """
rule gfa_linear_gfaffix:
    input:
        linear_rmac0_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize.rmac0.gfa"
    output:
        linear_gfaffix_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize.gfaffix.gfa",
        linear_gfaffix_info = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize.gfaffix.info"
    shell:
        """
        gfaffix -o {output.linear_gfaffix_gfa} \
            {input.linear_rmac0_gfa} \
            > {output.linear_gfaffix_info}
        """
        
rule remove_gfa_null_line:
    input:
        linear_gfaffix_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize.gfaffix.gfa"
    output:
        removed_null_linear_gfaffix_gfa = f"c7_graph_construction/chr_mc/{{chr}}/subgraph/subgraph{{subgraph_id}}/{config['prefix']}.{{chr}}.subgraph_{{subgraph_id}}.seqwish.smoothxg2.gfaffix.linearize.gfa"
    shell:
        """
        null_line=$(awk '{{if($0=="")print NR}}' {input.linear_gfaffix_gfa})
        awk -v null_line=$null_line '{{if(NR<null_line) print$0}}' {input.linear_gfaffix_gfa} > {output.removed_null_linear_gfaffix_gfa}
        """
