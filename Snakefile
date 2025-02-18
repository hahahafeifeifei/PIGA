configfile: "config/config.yaml"

# from functools import partial

with open('CHM13.20mb.interval') as f:
    intervals_list = [line.split("\t")[0].strip() for line in f]
config['intervals'] = intervals_list

# prefix = "t2t.grch38.58hifi.1064zmw"

rule all:
    input:
        expand("c9_diploid_path_infer/result/{sample}/{sample}.hap1.complete_assembly.polish.clip.fasta", sample=config['samples']),
        expand("c9_diploid_path_infer/result/{sample}/{sample}.hap2.complete_assembly.polish.clip.fasta", sample=config['samples'])
        
#load the included snakemake files.
# include: "rules/SR_var_calling.smk"
include: "rules/generate_personal_reference.smk"
include: "rules/split_minigraph.smk"
include: "rules/graph_construction.smk"
include: "rules/simplify_ml_pangenome.smk"
include: "rules/merge_pangenome.smk"
include: "rules/infer_diploid_path.smk"

#allow_missing=True