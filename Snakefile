import os

configfile: "config/config.yaml"
configfile: "config/tools.yaml"

TOOLS = {
    tool: path
     for tool, path in config['tools'].items()
}
# 
# # check the path of tools
# def validate_tool(tool_path):
#     
#     if not os.path.exists(tool_path):
#         raise ValueError(f"Tool path {tool_path} does not exist!")
#     return tool_path
# 
# # set a dictionary to check and store the path of tools.
# TOOLS = {
#     tool: validate_tool(path) for tool, path in config['tools'][category].items() for category in config['tools']
# }
# # from functools import partial

with open('CHM13.20mb.interval') as f:
    intervals_list = [line.split("\t")[0].strip() for line in f]
config['intervals'] = intervals_list

# prefix = "t2t.grch38.58hifi.1064zmw"

#load the included snakemake files.

include: "rules/utils.smk"
include: "rules/SR_var_calling.smk"
include: "rules/LR_var_calling.smk"
include: "rules/merge_snv.smk"
include: "rules/phase_snv.smk"
include: "rules/generate_personal_reference.smk"
include: "rules/draft_assembly.smk"
include: "rules/split_minigraph.smk"
include: "rules/graph_construction.smk"
include: "rules/simplify_ml_pangenome.smk"
include: "rules/merge_pangenome.smk"
include: "rules/infer_diploid_path.smk"

Call_SR_SNV = config.get("Call_SR_SNV", False)
Call_LR_SNV = config.get("Call_LR_SNV", False)
Merge_SNV = config.get("Merge_SNV", False)
Phase_SNV = config.get("Phase_SNV", False)
Generate_Personal_Reference = config.get("Generate_Personal_Reference", False)
Draft_Assembly = config.get("Draft_Assembly", False)
Split_Minigraph = config.get("Split_Minigraph", False)
Graph_Construction = config.get("Graph_Construction", False)
Simplify_ML_Pangenome = config.get("Simplify_ML_Pangenome", False)
Merge_Pangenome = config.get("Merge_Pangenome", False)
Infer_diploid_PATH = config.get("Infer_diploid_PATH", False)


if Call_SR_SNV:

    rule all:
        input:
            rules.all_SR_var_calling.input

elif Call_LR_SNV:

    rule all:
        input:
            rules.all_LR_var_calling.input
    
elif Merge_SNV:

    rule all:
        input:
            rules.all_merge_snv.input
    
elif Phase_SNV:

    rule all:
        input:
            rules.all_phase_snv.input
    
elif Generate_Personal_Reference:

    rule all:
        input:
            rules.all_generate_personal_reference.input
    
elif Draft_Assembly:

    rule all:
        input:
            rules.all_draft_assembly.input
    
elif Split_Minigraph:

    rule all:
        input:
            rules.all_split_minigraph.input
    
elif Graph_Construction:

    rule all:
        input:
            rules.all_graph_construction.input
    
elif Simplify_ML_Pangenome:
    
    rule all:
        input:
            rules.all_simplify_ml_pangenome.input
    
elif Merge_Pangenome:

    rule all:
        input:
            rules.all_merge_pangenome.input
    
elif Infer_diploid_PATH:

    rule all:
        input:
            rules.all_infer_diploid_path.input
else:

    rule all:
        input:
            rules.all_infer_diploid_path.input

        
#allow_missing=True
