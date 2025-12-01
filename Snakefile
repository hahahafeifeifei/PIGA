shell.prefix = ['set -e']

import os

configfile: "config/config.yaml"
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


samples_list, sex_list, sample_sex_dict = [], [], {}
with open(config['samples']) as f:
    for line in f:
        if not line.strip():
            continue

        sample = line.strip().split()[0]
        sex = line.strip().split()[1]

        samples_list.append(sample)
        sex_list.append(sex)
        sample_sex_dict[sample] = sex
config['sex'] = sample_sex_dict

wildcard_constraints:
    sample='|'.join(samples_list)


#load the included snakemake files.
include: "rules/utils.smk"


Call_SR_SNV = config.get("Call_SR_SNV", False)
Call_LR_SNV = config.get("Call_LR_SNV", False)
Merge_SNV = config.get("Merge_SNV", False)
Phase_SNV = config.get("Phase_SNV", False)
Generate_Personal_Reference = config.get("Generate_Personal_Reference", False)
Draft_Assembly = config.get("Draft_Assembly", False)
Consturct_Pangenome = config.get("Consturct_Pangenome", False)
Simplify_Pangenome = config.get("Simplify_Pangenome", False)
Merge_Pangenome = config.get("Merge_Pangenome", False)
Infer_Diploid_Path = config.get("Infer_Diploid_Path", False)


if Call_SR_SNV:

    include: "rules/call_sr_snv.smk"
    rule all:
        input:
            rules.all_call_sr_snv.input

elif Call_LR_SNV:

    include: "rules/call_lr_snv.smk"
    rule all:
        input:
            rules.all_call_lr_snv.input
    
elif Merge_SNV:

    include: "rules/merge_snv.smk"
    rule all:
        input:
            rules.all_merge_snv.input
    
elif Phase_SNV:

    include: "rules/phase_snv.smk"
    rule all:
        input:
            rules.all_phase_snv.input
    
elif Generate_Personal_Reference:

    include: "rules/generate_personal_reference.smk"
    rule all:
        input:
            rules.all_generate_personal_reference.input
    
elif Draft_Assembly:

    include: "rules/draft_assembly.smk"
    rule all:
        input:
            rules.all_draft_assembly.input
    
elif Consturct_Pangenome:

    include: "rules/construct_pangenome.smk"
    rule all:
        input:
            rules.all_construct_pangenome.input
    
elif Simplify_Pangenome:
    
    include: "rules/simplify_pangenome.smk"
    rule all:
        input:
            rules.all_simplify_pangenome.input
    
elif Merge_Pangenome:

    include: "rules/merge_pangenome.smk"
    rule all:
        input:
            rules.all_merge_pangenome.input
    
elif Infer_Diploid_Path:

    include: "rules/infer_diploid_path.smk"
    rule all:
        input:
            rules.all_infer_diploid_path.input
else:

    include: "rules/call_sr_snv.smk"
    include: "rules/call_lr_snv.smk"
    include: "rules/merge_snv.smk"
    include: "rules/phase_snv.smk"
    include: "rules/generate_personal_reference.smk"
    include: "rules/draft_assembly.smk"
    include: "rules/construct_pangenome.smk"
    include: "rules/simplify_pangenome.smk"
    include: "rules/merge_pangenome.smk"
    include: "rules/infer_diploid_path.smk"
    rule all:
        input:
            rules.all_infer_diploid_path.input

        
#allow_missing=True
