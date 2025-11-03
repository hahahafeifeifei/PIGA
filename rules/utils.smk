def get_sr_input_fastqs(wildcards):
    return config["sr_fastqs"]

def get_zmw_input_fastqs(wildcards):
    return config["lr_zmw_fastqs"]

def get_hifi_input_fastqs(wildcards):
    return config["lr_hifi_fastqs"]

def get_sex(wildcards):
    sex = config['sex'][wildcards.sample]
    return sex

def get_sex_specific_chr_list(wildcards):
    sex = config['sex'][wildcards.sample]
    
    chr_list = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    if sex == 'male':
        return chr_list + ['chrY']
    return chr_list

def get_sr_vcf_input(wildcards):
    if "sr_vcf" in config:
        return config["sr_vcf"]
    else:
        return f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.vcf.gz"

def get_lr_vcf_input(wildcards):
    if "lr_vcf" in config:
        return config["lr_vcf"]
    else:
        return f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.beagle.vcf.gz"

def get_merge_merfin_filter_vcf_input(wildcards):
    if "merge_merfin_filter_vcf" in config:
        return config["merge_merfin_filter_vcf"]
    else:
        return f"c3_merge_snv/merged_vcf/{config['prefix']}.consensus.merfin.vcf.gz"

def get_zmw_bam_input(wildcards):
    if "zmw_bam" in config:
        return config["zmw_bam"]
    else:
        return "c2_call_lr_snv/sample_bam/{sample}/{sample}.zmw.srt.bam"

def get_consensus_fasta_input(wildcards):
    if "consensus_fasta" in config:
        return config["consensus_fasta"]
    else:
        return "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.fasta"

def get_chain_input(wildcards):
    if "chain" in config:
        return config["chain"]
    else:
        return "c5_personal_ref/sample_reference/{sample}/{sample}.personal_ref.chain"

def get_sample_vcf_input(wildcards):
    if "sample_vcf" in config:
        return config["sample_vcf"]
    else:
        return "c4_phase_snv/sample_vcf/{sample}/{sample}.shapeit.vcf.gz"

def get_sample_meryl_input(wildcards):
    if "sample_meryl" in config:
        return config["sample_meryl"]
    else:
        return "c3_merge_snv/sample_meryl/{sample}/{sample}.meryl"

def concat_final_phase_vcf_sex_specific_chrlist(wildcards):
    sex = config['sex'][wildcards.sample]
    chr_list = [f'chr{i}' for i in range(1, 23)]
    if sex == 'female':
        return chr_list + ['chrX']
    return chr_list


#for a given chromosome, generate all subgraph_id.
def get_all_subgraph_assembly_gfa_files(wildcards, prefix):
    if "subgraph_list" in config:
        chr_subgraph_combination_file = config["subgraph_list"]
    else:
        chr_subgraph_combination = checkpoints.prepare_subgraph_list.get().output[0]
    chr_pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, subgraph_id = line.strip().split("\t")
            if chr == wildcards.chr:
                if "subgraph_assembly_gfa" in config:
                    subgraph_assembly_gfa_file = config['subgraph_assembly_gfa'].format(chr = wildcards.chr, subgraph_id = wildcards.subgraph_id)
                else:
                    subgraph_assembly_gfa_file = f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.assembly.gfa"
                
                chr_pairs.append(subgraph_assembly_gfa_file)
    return chr_pairs

def get_all_subgraph_variant_path_files(wildcards, prefix):
    if "subgraph_list" in config:
        chr_subgraph_combination_file = config["subgraph_list"]
    else:
        chr_subgraph_combination = checkpoints.prepare_subgraph_list.get().output[0]
    chr_pairs = []
    with open(chr_subgraph_combination) as f:
        for line in f:
            chr, subgraph_id = line.strip().split("\t")
            if chr == wildcards.chr:
                if "subgraph_variant_path" in config:
                    subgraph_variant_path_file = config['subgraph_variant_path'].format(chr = wildcards.chr, subgraph_id = wildcards.subgraph_id)
                else:
                    subgraph_variant_path_file = f"c7_graph_construction/chr_mc/{chr}/subgraph/subgraph{subgraph_id}/{prefix}.{chr}.subgraph_{subgraph_id}.seqwish.smoothxg2.gfaffix.linearize.TVR90.variant_project.gfaffix.chop.ids.variant.path"
                
                chr_pairs.append(subgraph_assembly_gfa_file)
    return chr_pairs

def get_hapl_input(wildcards):
    if "hapl" in config:
        return config["hapl"]
    else:
        return "c7_graph_construction/graph_merge/CKCG.merge.assembly.hapl"

def get_gbz_input(wildcards):
    if "gbz" in config:
        return config["gbz"]
    else:
        return "c7_graph_construction/graph_merge/CKCG.merge.assembly.gbz"

def get_gfa_input(wildcards):
    if "gfa" in config:
        return config["gfa"]
    else:
        return "c7_graph_construction/graph_merge/CKCG.{chr}.assembly.gfa"

def get_variant_path_input(wildcards):
    if "variant_path" in config:
        return config["variant_path"]
    else:
        return "c7_graph_construction/graph_merge/CKCG.{chr}.variant.path"
        
def get_hap1_adaptor_masked_fa_input(wildcards):
    if "hap1_adaptor_masked_fa" in config:
        return config["hap1_adaptor_masked_fa"]
    else:
        return "c6_draft_assembly/result/{sample}/assembly/{sample}.hap1.adaptor_masked.fasta"

def get_hap2_adaptor_masked_fa_input(wildcards):
    if "hap2_adaptor_masked_fa" in config:
        return config["hap2_adaptor_masked_fa"]
    else:
        return "c6_draft_assembly/result/{sample}/assembly/{sample}.hap2.adaptor_masked.fasta"
