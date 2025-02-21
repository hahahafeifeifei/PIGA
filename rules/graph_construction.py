chr_list = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


#prefix = '1kcp'
            
def get_linear_gfaffix_gfa_files(wildcards, prefix):
    chr_subgraph_combination = checkpoints.check_chr_subgraph_combination.get().output[0]
    pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr_id, subgraph_id = line.strip().split("\t")
            pairs.append(f"c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfa")
    return pairs
            
# rule all:
#     input:
#         partial(get_linear_gfaffix_gfa_files, prefix = config['prefix'])
    

#TODO: the checked list is not necessary.Maybe edit the logic in the future.
checkpoint check_chr_subgraph_combination:
    input:
        subgraph_list = "c8_graph_construction/chr_mc/subgraph.list" 
    output:
        list_checked = "c8_graph_construction/chr_mc/subgraph.list.checked" 
    run:
        with open(input[0], "r") as f_in, open(output[0], "w") as f_out:
            for line in f_in:
                chr_id, subgraph_id = line.strip().split()
                f_out.write(f"{chr_id}\t{subgraph_id}\n")    



#Here subgraph_seqfile only to show the graph splitting has been finished.
rule prepare_subgraph_list:
    input:
        subgraph_seqfile = expand("c8_graph_construction/chr_mc/{chr}/subgraph/subgraph0/{prefix}.{chr}.subgraph_0.seqfile", chr=chr_list, prefix=config['prefix'])
    output:
        subgraph_list = "c8_graph_construction/chr_mc/subgraph.list" 
    shell:
        """
        for chr in `ls -d chr*`;do  for subgraph_id in `ls -d ${{chr}}/subgraph/subgraph*`;do idx=$(echo ${{subgraph_id}} | awk '{{split($1,a,"/");print a[length(a)]}}' | awk -v FS='subgraph' '{{print $2}}'); echo -e ${{chr}}"\\t"${{idx}}; done; done > c8_graph_construction/chr_mc/subgraph.list
        """
    
    
    
rule subgraph_seqwish:
    input:
        paf = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.paf"
        
    output:
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.gfa"
    threads: 1
    params:
        prefix = config['prefix'],
        chr_subgraph_dir = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}"
    shell:
        """
        cat {params.chr_subgraph_dir}/fa/*.fasta > {params.chr_subgraph_dir}/{params.prefix}.{chr}.subgraph_{wildcards.subgraph_id}.fasta
        seqwish -P -t {threads} \
        -s {params.chr_subgraph_dir}/{params.prefix}.{chr}.subgraph_{wildcards.subgraph_id}.fasta \
        -p {params.chr_subgraph_dir}/{params.prefix}.{chr}.subgraph_{wildcards.subgraph_id}.paf \
        -g {params.chr_subgraph_dir}/{params.prefix}.{chr}.subgraph_{wildcards.subgraph_id}.seqwish.origin.gfa
        
        awk -v OFS='\t' '{{if(substr($1,1,1)=="P") {{split($2,a,"id=");split(a[2],b,"|");if(b[1]=="CHM13" || b[1]=="GRCh38" || b[1]=="_MINIGRAPH_") name=b[1]"#"b[2];else {{split(b[1],c,".");name=c[1]"#"c[2]-1"#"b[2]"#0"}}  print $1,name,$3,$4 }}else  print$0  }}' {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.origin.gfa > {params.chr_subgraph_dir}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.gfa
        
        rm {params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.fasta
        rm {params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.origin.gfa
        """

        
rule subgraph_smoothxg_1:
    input:
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.gfa"
    output:
        smoothxg_1_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg1.gfa"
    threads: 1
    params:
        chr_subgraph_dir = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}"
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
        smoothxg_1_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg1.gfa"
        
    output:
        smoothxg_2_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfa"
    threads: 1
    params:
        chr_subgraph_dir = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}"
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
#         smoothxg_2_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfa"
        
#     output:
#         smoothxg_3_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg3.gfa"
#     params:
#         chr_subgraph_dir = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}"
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
        smoothxg_2_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfa" 
    output:
        gfaffix_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.gfa",
        gfaffix_info = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.info"
    shell:
        """
        gfaffix \
        {input.smoothxg_2_gfa} \
        -o {output.gfaffix_gfa} \
        > {output.gfaffix_info}
        
        """

rule gfa_format_norm:
    input:
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.gfa"    
    output:
        norm_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.norm.gfa"
    params:
        prefix = config['prefix']
    shell:
        """
        vg convert -gf {input.gfa} > c8_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.W.gfa
        
        awk -v OFS='\t' '{{ if(substr($0,1,1)=="W" && $2!="_MINIGRAPH_") {{split($4,a,":");split(a[2],b,"-");print$1,$2,$3,a[1],$5+b[1]-1,$6+b[1]-1,$7}} else print$0 }}' \
        c8_graph_construction/chr_mc/{wildcards.chr}/subgraph/subgraph{wildcards.subgraph_id}/{params.prefix}.{wildcards.chr}.subgraph_{wildcards.subgraph_id}.seqwish.smoothxg2.gfaffix.W.gfa \
        > {output.norm_gfa}
        
        """
        
rule record_gfa_cycle_region:
    input:
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.norm.gfa"
    output:
        cycle_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.norm.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/cycle_region.py \
        {input.gfa} \
        {output.cycle_region}
        """
        

rule record_gfa_high_cov_region:
    input:
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.norm.gfa"
    output:
        high_cov_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.high_coverage.bed"
    shell:
        """
        python3 scripts/graph-simplification/gfa_high-coverage_bed.py \
        {input.gfa} \
        45200 \
        200 \
        {output.high_cov_region}
        """        
        
#TODO: how to generate the mask region? 
rule gfa_regions_merge:
    input:
        mask="/storage/yangjianLab/wangyifei/project/01.CKCG/07.CLR_Pangenome/graph_construction/cactus/clip_bed/mask_bed/{chr}.mask.bed",
        high_cov_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.high_coverage.bed"
    output:
        merged_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.merge.bed"
    shell:
        """
        cat {mask} {input.high_cov_region} | bedtools sort | bedtools merge -d 1000 > {output.merged_region}
        
        """        
rule gfa_filter:
    input:
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.norm.gfa",
        merged_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.merge.bed"
    output:
        filtered_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.clip.gfa" 
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
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.clip.gfa"
    output:
        cycle_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.clip.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/cycle_region.py \
        {input.gfa} \
        {output.cycle_region}
        """

#TODO: the include and the candidate?
rule gfa_linear_1:
    input:
        gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.clip.gfa"
    output:
        linear_1_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize1.gfa",
        linear_1_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize1.region",
        linear_1_cycle_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize1.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub_muti_match.py \
        {input.gfa} 5000 50 ${include} ${candidate} \
        {output.linear_1_gfa} \
        {output.linear_1_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
        {output.linear_1_gfa} \
        {output.linear_1_cycle_region}
        """
        
rule gfa_linear_2:
    input:
        linear_1_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize1.gfa"
    output:
        linear_2_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize2.gfa",
        linear_2_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize2.region",
        linear_2_cycle_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize2.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub.py \
        {input.linear_1_gfa} \
        5000 50 10000 ${include} ${candidate} \
        {output.linear_2_gfa} \
        {output.linear_2_cycle_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
        {output.linear_2_gfa} \
        {output.linear_2_cycle_region}
        """
        
rule gfa_linear_3:
    input:
        linear_2_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize2.gfa"
    output:
        linear_3_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize3.gfa",
        linear_3_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize3.region",
        linear_3_cycle_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize3.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub.py \
        {input.linear_2_gfa} \
        7500 50 5000 ${include} ${candidate} \
        {output.linear_3_gfa} \
        {output.linear_3_cycle_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
        {output.linear_3_gfa} \
        {output.linear_3_cycle_region}
        """
        
rule gfa_linear_4:
    input:
        linear_3_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize3.gfa"
    output:
        linear_4_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize4.gfa",
        linear_4_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize4.region",
        linear_4_cycle_region = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize4.cycle.region"
    shell:
        """
        python3 scripts/graph-simplification/gfa_untangle_unsub.py \
        {input.linear_3_gfa} \
        9000 50 2000 ${include} ${candidate} \
        {output.linear_4_gfa} \
        {output.linear_4_cycle_region}
        
        python3 scripts/graph-simplification/cycle_region.py \
        {output.linear_4_gfa} \
        {output.linear_4_cycle_region}
        """

        
        
        
rule gfa_linear_rmac0:
    input:
        linear_4_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize4.gfa"
    output:
        linear_rmac0_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.rmac0.gfa"
    shell:
        """
        python3 scripts/graph-simplification/gfa_remove_ac0.py \
        {input.linear_4_gfa} \
        {output.linear_rmac0_gfa}
        """
rule gfa_linear_gfaffix:
    input:
        linear_rmac0_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.rmac0.gfa"
    output:
        linear_gfaffix_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfaffix.gfa",
        linear_gfaffix_info = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfaffix.info"
    shell:
        """
        gfaffix -o {output.linear_gfaffix_gfa} \
        {input.linear_rmac0_gfa} \
        > {output.linear_gfaffix_info}
        """
        
rule remove_gfa_null_line:
    input:
        linear_gfaffix_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfaffix.gfa"
    output:
        removed_null_linear_gfaffix_gfa = "c8_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{config['prefix']}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.gfa"
    shell:
        """
        null_line=$(awk '{{if($0=="")print NR}}' {input.linear_gfaffix_gfa})
        awk -v null_line=$null_line '{{if(NR<null_line) print$0}}' {input.linear_gfaffix_gfa} > {output.removed_null_linear_gfaffix_gfa}
        """