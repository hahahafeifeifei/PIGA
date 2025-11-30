chr_list= [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

def get_all_subgraph_files(wildcards):
    subgraph_id_file = checkpoints.minigraph_aln_partition.get().output.subgraph_id_list
    with open(subgraph_id_file) as f:
        subgraph_ids = [line.strip() for line in f if line.strip()]
    return expand(f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.gfa", id=subgraph_ids)

rule all_construct_pangenome:
    input:
        get_all_subgraph_files

rule rename_reference:
    input:
        chm13_fa=config["reference"]["CHM13"],
        grch38_fa=config["reference"]["GRCh38"]
    output:
        chm13_rename_fa='c7_graph_construction/fasta/CHM13.fasta',
        chm13_rename_fai='c7_graph_construction/fasta/CHM13.fasta.fai',
        grch38_rename_fa='c7_graph_construction/fasta/GRCh38.fasta',
        grch38_rename_fa_fai='c7_graph_construction/fasta/GRCh38.fasta.fai'
    resources:
        max_mem_gb=30,
    shell:
        """
        awk '{{if(substr($0,1,1)==">") print ">CHM13."substr($1,2,length($1));else print $0}}' {input.chm13_fa} > {output.chm13_rename_fa}
        samtools faidx {output.chm13_rename_fa}
        awk '{{if(substr($0,1,1)==">") print ">GRCh38."substr($1,2,length($1));else print $0}}' {input.grch38_fa} > {output.grch38_rename_fa}
        samtools faidx {output.grch38_rename_fa}
        """

external_assembly_id_dict, external_assembly_id_list = {}, []
with open(get_external_assembly_list()) as f:
    for line in f:
        if not line.strip():
            continue
        assembly_id = line.strip().split()[0]
        assembly = line.strip().split()[1]
        external_assembly_id_list.append(assembly_id)
        external_assembly_id_dict[assembly_id] = assembly

rule rename_external_assembly:
    input:
        fa = get_external_assembly_fa
    output:
        rename_fa = 'c7_graph_construction/fasta/{external_assembly_id}.fasta',
        rename_fai = 'c7_graph_construction/fasta/{external_assembly_id}.fasta.fai'
    wildcard_constraints:
        external_assembly_id = "|".join([re.escape(str(i)) for i in external_assembly_id_list])
    resources:
        max_mem_gb=30
    shell:
        """
        sample=$(echo {wildcards.external_assembly_id} | awk '{{split($1,a,".");print a[1]}}')
        haplotype=$(echo {wildcards.external_assembly_id} | awk '{{split($1,a,".");print a[2]}}')
        awk -v sample=$sample -v hap=$haplotype '{{if(substr($0,1,1)==">") print ">"sample"."hap"."substr($1,2,length($1));else print $0}}' {input.fa} > {output.rename_fa}
        samtools faidx {output.rename_fa}
        """

rule external_minigraph:
    input:
        chm13_fa='c7_graph_construction/fasta/CHM13.fasta',
        grch38_fa='c7_graph_construction/fasta/GRCh38.fasta',
        external_assembly_fa = expand('c7_graph_construction/fasta/{external_assembly_id}.fasta', external_assembly_id = external_assembly_id_list)
    output:
        minigraph_external_gfa = f"c7_graph_construction/{config['prefix']}.minigraph_external.gfa"
    threads: 16
    resources:
        max_mem_gb=300
    shell:
        """
        minigraph -c -x ggs -l 10000 -t {threads} {input.chm13_fa} {input.grch38_fa} {input.external_assembly_fa} > {output.minigraph_external_gfa}
        """

rule external_chr_minigraph:
    input:
        minigraph_external_gfa = f"c7_graph_construction/{config['prefix']}.minigraph_external.gfa"
    output:
        minigraph_external_chr_gfa = f"c7_graph_construction/chr_gfa/{config['prefix']}.minigraph_external.{{chr}}.gfa"
    wildcard_constraints:
        chr = "|".join(chr_list)
    threads: 4
    resources:
        max_mem_gb=30
    shell:
        """
        if [ {wildcards.chr} == "chrM" ]
        then
            grep chrM {input.minigraph_external_gfa} > {output.minigraph_external_chr_gfa}
        else
            gfatools view -R CHM13.{wildcards.chr} {input.minigraph_external_gfa} > {output.minigraph_external_chr_gfa}
        fi
        """

internal_assembly_id_dict, internal_assembly_id_list = {}, []
with open(get_internal_assembly_list()) as f:
    for line in f:
        if not line.strip():
            continue
        assembly_id = line.strip().split()[0]
        assembly = line.strip().split()[1]
        internal_assembly_id_list.append(assembly_id)
        internal_assembly_id_dict[assembly_id] = assembly

rule rename_internal_assembly:
    input:
        fa = get_internal_assembly_fa 
    output:
        rename_fa = 'c7_graph_construction/fasta/{internal_assembly_id}.fasta',
        rename_fai = 'c7_graph_construction/fasta/{internal_assembly_id}.fasta.fai'
    wildcard_constraints:
        internal_assembly_id = "|".join([re.escape(str(i)) for i in internal_assembly_id_list])
    resources:
        max_mem_gb=30
    shell:
        """
        sample=$(echo {wildcards.internal_assembly_id} | awk '{{split($1,a,".");print a[1]}}')
        haplotype=$(echo {wildcards.internal_assembly_id} | awk '{{split($1,a,".");print a[2]+1}}')
        awk -v sample=$sample -v hap=$haplotype '{{if(substr($0,1,1)==">") print ">"sample"."hap"."substr($1,2,length($1));else print $0}}' {input.fa} > {output.rename_fa}
        samtools faidx {output.rename_fa}
        """

rule split_internal_assembly_by_chr:
    input:
        fa = 'c7_graph_construction/fasta/{internal_assembly_id}.fasta'
    output:
        chr_fa = "c7_graph_construction/chr_fasta/{chr}/{internal_assembly_id}.{chr}.fasta"
    wildcard_constraints:
        chr = "|".join(chr_list)
    threads: 1
    resources:
        max_mem_gb = 30
    shell:
        """
        seqkit grep -w 0 -r -p {wildcards.chr}- {input.fa} > {output.chr_fa}
        """

rule internal_chr_minigraph:
    input:
        minigraph_external_chr_gfa = f"c7_graph_construction/chr_gfa/{config['prefix']}.minigraph_external.{{chr}}.gfa",
        internal_assembly_fa = expand('c7_graph_construction/chr_fasta/{{chr}}/{internal_assembly_id}.{{chr}}.fasta', internal_assembly_id = internal_assembly_id_list)
    output:
        minigraph_chr_gfa = f"c7_graph_construction/chr_gfa/{config['prefix']}.minigraph.{{chr}}.gfa"
    threads: 8
    resources:
        max_mem_gb = 60
    shell:
        """
        minigraph -c -x ggs -l 10000 -t {threads} {input.minigraph_external_chr_gfa} {input.internal_assembly_fa} | sed "s/s//g" | vg convert -fg - > {output.minigraph_chr_gfa}
        """

rule concat_chr_minigraph:
    input:
        minigraph_external_chr_gfas = expand(f"c7_graph_construction/chr_gfa/{config['prefix']}.minigraph.{{chr}}.gfa", chr=chr_list)
    output:
        minigraph_gfa = f"c7_graph_construction/{config['prefix']}.minigraph.gfa"
    threads: 16
    resources:
        max_mem_gb = 200
    shell:
        """
        vg concat {input.minigraph_external_chr_gfas} > {output.minigraph_gfa}
        """

rule minigraph_clip:
    input:
        minigraph_gfa = f"c7_graph_construction/{config['prefix']}.minigraph.gfa"
    output:
        snarls=f"c7_graph_construction/{config['prefix']}.minigraph.snarls",
        node_clip_gfa=f"c7_graph_construction/{config['prefix']}.minigraph.node_clip.gfa",
        edge_clip_gfa=f"c7_graph_construction/{config['prefix']}.minigraph.node_edge_clip.gfa",
        node_clip_rgfa=f"c7_graph_construction/{config['prefix']}.minigraph.node_clip.rgfa",
        node_len=f"c7_graph_construction/{config['prefix']}.minigraph.node_len.tsv",
        subgraph_info=f"c7_graph_construction/{config['prefix']}.subgraph.info"
    params:
        tmp_dir = f"c7_graph_construction"
    threads: 16
    resources:
        max_mem_gb = 200
    shell:
        """
        vg snarls {input.minigraph_gfa} | vg view -R - > {output.snarls}
        python3 scripts/construct_pangenome/gfa_border_node_select.py {input.minigraph_gfa} {output.snarls} {params.tmp_dir} {output.node_clip_gfa} {output.edge_clip_gfa} {output.subgraph_info}
        vg convert -g {output.node_clip_gfa} -f -Q CHM13.chr | awk '{{if(substr($0,1,1)=="S" && $6!="SR:i:0")print$0"\\tSN:Z:Other\\tSO:i:0\\tSR:i:1";else print$0}}' | \
        awk -v OFS='\\t' '{{if(substr($0,1,1)=="S")print$1,"s"$2,$3,$4,$5,$6;else {{if(substr($0,1,1)=="L")print$1,"s"$2,$3,"s"$4,$5,$6;else print$0}} }}' > {output.node_clip_rgfa}
        grep ^S {output.node_clip_rgfa} | awk '{{print $2"\\t"length($3)}}' > {output.node_len}
        """


#Rule: minigraph alignment
rule minigraph_alignment_reference:
    input:
        node_clip_rgfa = f"c7_graph_construction/{config['prefix']}.minigraph.node_clip.rgfa",
        chm13_fa = 'c7_graph_construction/fasta/CHM13.fasta',
        grch38_fa = 'c7_graph_construction/fasta/GRCh38.fasta'
    output:
        chm13_gaf = "c7_graph_construction/gaf/CHM13.gaf",
        grch38_gaf = "c7_graph_construction/gaf/GRCh38.gaf"
    resources:
        max_mem_gb = 100
    threads: 10
    shell:
        """
        minigraph -x asm -t {threads} --vc -c {input.node_clip_rgfa} {input.chm13_fa} | gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > {output.chm13_gaf}
        minigraph -x asm -t {threads} --vc -c {input.node_clip_rgfa} {input.grch38_fa} | gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > {output.grch38_gaf}
        """

rule minigraph_alignment_external_assembly:
    input:
        node_clip_rgfa = f"c7_graph_construction/{config['prefix']}.minigraph.node_clip.rgfa",
        fa = 'c7_graph_construction/fasta/{external_assembly_id}.fasta'
    output:
        gaf = "c7_graph_construction/gaf/{external_assembly_id}.gaf"
    wildcard_constraints:
        external_assembly_id = "|".join([re.escape(str(i)) for i in external_assembly_id_list])
    resources:
        max_mem_gb = 100
    threads: 10
    shell:
        """
        minigraph -x asm -t {threads} --vc -c {input.node_clip_rgfa} {input.fa} | gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > {output.gaf}
        """

rule minigraph_alignment_internal_assembly:
    input:
        node_clip_rgfa = f"c7_graph_construction/{config['prefix']}.minigraph.node_clip.rgfa",
        fa = 'c7_graph_construction/fasta/{internal_assembly_id}.fasta'
    output:
        gaf = "c7_graph_construction/gaf/{internal_assembly_id}.gaf"
    wildcard_constraints:
        internal_assembly_id = "|".join([re.escape(str(i)) for i in internal_assembly_id_list])
    resources:
        max_mem_gb = 100
    threads: 10
    shell:
        """
        minigraph -x asm -t {threads} --vc -c {input.node_clip_rgfa} {input.fa} | gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > {output.gaf}
        """

checkpoint minigraph_aln_partition:
    input:
        gafs = expand("c7_graph_construction/gaf/{assembly_id}.gaf", assembly_id = ["GRCh38", "CHM13"] + internal_assembly_id_list + external_assembly_id_list),
        node_edge_clip_gfa = f"c7_graph_construction/{config['prefix']}.minigraph.node_edge_clip.gfa",
        node_len = f"c7_graph_construction/{config['prefix']}.minigraph.node_len.tsv"
    output:
        subgraph_dir = directory("c7_graph_construction/subgraph"),
        gaf = f"c7_graph_construction/{config['prefix']}.gaf",
        paf = f"c7_graph_construction/{config['prefix']}.paf",
        subgraph_id_list = f"c7_graph_construction/subgraph_id.list",
    params:
        prefix = f"c7_graph_construction/subgraph/{config['prefix']}_subgraph"
    resources:
        max_mem_gb=100,
        runtime_hrs=30
    threads: 8
    shell:
        """
        mkdir {output.subgraph_dir}
        cat {input.gafs} > {output.gaf}
        gaf2paf -l {input.node_len} {output.gaf} > {output.paf}
        vg chunk -C -x {input.node_edge_clip_gfa} --prefix {params.prefix} -O gfa
        python3 scripts/construct_pangenome/subgraph_paf_fa.py {params.prefix} {output.paf} {output.gaf} {output.subgraph_dir}
        ls {params.prefix}_*.gfa | awk '{{split($1,a,"_");split(a[length(a)],b,".");print b[1]}}' > {output.subgraph_id_list}
        """

#Rule: minigraph sequence partition
rule minigraph_seq_partition:
    input:
        fa_dir = "c7_graph_construction/fasta",
        bed = f"c7_graph_construction/subgraph/{config['prefix']}_subgraph_{{id}}.bed",
        gfa = f"c7_graph_construction/subgraph/{config['prefix']}_subgraph_{{id}}.gfa", 
        paf = f"c7_graph_construction/subgraph/{config['prefix']}_subgraph_{{id}}.paf"
    output:
        subgraph_dir=directory("c7_graph_construction/subgraph/subgraph_{id}"),
        subgraph_fa_dir=directory("c7_graph_construction/subgraph/subgraph_{id}/fasta"),
        fa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.fasta",
        bed = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.bed",
        paf = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.paf",
        gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.gfa"
    params:
        prefix = f"c7_graph_construction/subgraph/{config['prefix']}_subgraph"
    wildcard_constraints:
        id='[0-9]+'
    resources:
        max_mem_gb=100,
        runtime_hrs=30
    threads: 1
    shell:
        """
        mkdir -p {output.subgraph_dir}
        mkdir -p {output.subgraph_fa_dir}
        awk '{{print $1}}' {input.bed} | while read sample
        do
            awk -v sample=$sample '{{if($1==sample)print$2}}' {input.bed} > {output.subgraph_fa_dir}/$sample.subgraph_{wildcards.id}.bed
            samtools faidx -r {output.subgraph_fa_dir}/$sample.subgraph_{wildcards.id}.bed {input.fa_dir}/$sample.fasta > {output.subgraph_fa_dir}/$sample.subgraph_{wildcards.id}.fasta
        done
        awk '{{if($1=="S")print">_MINIGRAPH_.s"$2"\\n"$3}}' {input.gfa} > {output.subgraph_fa_dir}/_MINIGRAPH_.subgraph_{wildcards.id}.fasta
        cat {output.subgraph_fa_dir}/*fasta > {output.fa}
        mv {input.bed} {output.bed}
        mv {input.gfa} {output.gfa}
        awk -v OFS='\\t' '{{$6="_MINIGRAPH_."$6; print $0}}' {input.paf} > {output.paf}
        rm {input.paf}
        """


rule seqwish:
    input:
        fa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.fasta",
        paf = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.paf"
    output:
        raw_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.origin.gfa",
        seqwish_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.gfa"
    params:
        subgraph_dir = 'c7_graph_construction/subgraph/subgraph_{id}'
    resources:
        max_mem_gb=lambda wildcards, attempt: 100 * attempt,
    threads: 8
    shell:
        """
        seqwish -P -t {threads} -s {input.fa} -p {input.paf} -g {output.raw_gfa} -b {params.subgraph_dir}
        awk -v OFS='\\t' '{{if(substr($1,1,1)=="P") {{split($2,a,".");if(a[1]=="CHM13" || a[1]=="GRCh38" || a[1]=="_MINIGRAPH_") name=a[1]"#"a[1]"."a[2];else {{name=a[1]"#"a[2]"#"a[3]"#0"}}  print $1,name,$3,$4 }}else  print$0 }}' {output.raw_gfa} > {output.seqwish_gfa}
        """


rule smoothxg:
    input:
        fa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.fasta",
        seqwish_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.gfa"
    output:
        smoothxg_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfa"
    params:
        smoothxg_dir = 'c7_graph_construction/subgraph/subgraph_{id}/tmp'
    resources:
        max_mem_gb=lambda wildcards, attempt: 100 * attempt,
        runtime_hrs=lambda wildcards, attempt: 20 * attempt
    threads: 8
    shell:
        """
        mkdir -p {params.smoothxg_dir}
        sample_number=$(grep ">" {input.fa} | awk '{{split($1,a,".");if(length(a)==2) print a[1];else print a[1]"."a[2] }}' | sort -u | wc -l)
        smoothxg -t {threads} -g {input.seqwish_gfa} --base {params.smoothxg_dir} --chop-to 100 -I 0.98 -R 0 -j 0 -e 0 -l 1400,2200 -p 1,4,6,2,26,1 -O 0.001 -Y $[100*$sample_number] -d 0 -D 0 -V -c 200M -W 1 -o {output.smoothxg_gfa}
        """

        
rule gfaffix:
    input:
        smoothxg_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfa"
    output:
        gfaffix_gfa = f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.gfa"
    resources:
        max_mem_gb=lambda wildcards, attempt: 100 * attempt,
        runtime_hrs=lambda wildcards, attempt: 20 * attempt
    threads: 8
    shell:
        """
        gfaffix {input.smoothxg_gfa} | vg convert -fg - | \
        awk -v OFS='\\t' '{{if(substr($0,1,1)=="W" && $2!="_MINIGRAPH_") {{split($4,a,":");split(a[2],b,"-");print$1,$2,$3,a[1],$5+b[1]-1,$6+b[1]-1,$7}} else print$0 }}' > {output.gfaffix_gfa}
        """


