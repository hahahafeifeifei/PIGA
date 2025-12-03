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

def get_merge_phased_vcf_input(wildcards):
    if "phased_snv" in config:
        return config["phased_snv"]
    else:
        return f"c4_phase_snv/merged_vcf/{config['prefix']}.consensus.whatshap.shapeit.vcf.gz"

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

def get_sample_phased_vcf_input(wildcards):
    if "sample_phased_vcf" in config:
        return config["sample_phased_vcf"]
    else:
        return "c4_phase_snv/sample_vcf/{sample}/{sample}.shapeit.vcf.gz"

def get_sample_meryl_input(wildcards):
    if "sample_meryl" in config:
        return config["sample_meryl"]
    else:
        return "c3_merge_snv/sample_meryl/{sample}/{sample}.meryl/merylIndex"

def concat_final_phase_vcf_sex_specific_chrlist(wildcards):
    sex = config['sex'][wildcards.sample]
    chr_list = [f'chr{i}' for i in range(1, 23)]
    if sex == 'female':
        return chr_list + ['chrX']
    return chr_list


def get_external_assembly_fa(wildcards):
    return external_assembly_id_dict[wildcards.external_assembly_id]


def generate_internal_assembly_dict(config):
    
    if "internal_assembly_list" in config:
        internal_assembly_list_file = config["internal_assembly_list"]
    else:
        ### Step by step mode, 'c6_draft_assembly/sample_assembly/internal_assembly.list' has been provided.
        if os.path.exists('c6_draft_assembly/sample_assembly/internal_assembly.list'):
            internal_assembly_list_file = 'c6_draft_assembly/sample_assembly/internal_assembly.list'
        ### Run All mode, 'c6_draft_assembly/sample_assembly/internal_assembly.list' hasn't been generated yet.
        else:
            internal_assembly_list_file = checkpoints.assembly_list.get().output.internal_assembly_list


    with open(internal_assembly_list_file) as f:

      internal_assembly_dict = {}
      for line in f:
          if not line.strip():
              continue
          assembly_id = line.strip().split()[0]
          assembly = line.strip().split()[1]
          internal_assembly_dict[assembly_id] = assembly
      return internal_assembly_dict

def get_internal_assembly_id_list(wildcards):
    
    if "internal_assembly_dict" not in config:
        config['internal_assembly_dict'] = generate_internal_assembly_dict(config)

    return [id for id in config['internal_assembly_dict'].keys()]

def get_internal_assembly_fa(wildcards):
    
    if "internal_assembly_dict" not in config:
        config['internal_assembly_dict'] = generate_internal_assembly_dict(config)

    return config['internal_assembly_dict'][wildcards.internal_assembly_id]

def get_already_subgraph_ids(wildcards):
    if "subgraph_id_list" in config:
        subgraph_id_list = config["subgraph_id_list"]
    else:
        ### Step by step mode, 'c7_graph_construction/subgraph_id_list' has been provided.
        if os.path.exists('c7_graph_construction/subgraph_id.list'):
            subgraph_id_list = 'c7_graph_construction/subgraph_id.list'
        ### Run All mode, 'c7_graph_construction/subgraph_id_list' hasn't been generated yet.
        else:
            subgraph_id_list = checkpoints.minigraph_aln_partition.get().output.subgraph_id_list
    with open(subgraph_id_list) as f:
        return [line.strip() for line in f if line.strip()]

def get_subgraph_fa(wildcards):
    if "subgraph_fa" in config:
        return config["subgraph_fa"]
    else:
        return f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.fasta"

def get_subgraph_gfa(wildcards):
    if "subgraph_gfa" in config:
        return config["subgraph_gfa"]
    else:
        return f"c7_graph_construction/subgraph/subgraph_{{id}}/{config['prefix']}_subgraph_{{id}}.seqwish.smoothxg.gfaffix.gfa"

def get_hapl_input(wildcards):
    if "hapl" in config:
        return config["hapl"]
    else:
        return f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.hapl"

def get_gbz_input(wildcards):
    if "gbz" in config:
        return config["gbz"]
    else:
        return f"c7_graph_construction/graph_merge/{config['prefix']}.merge.assembly.gbz"

def get_gfa_input(wildcards):
    if "gfa" in config:
        return config["gfa"]
    else:
        return f"c7_graph_construction/graph_merge/{config['prefix']}.nopath.gfa"

def get_variant_path_input(wildcards):
    if "variant_path" in config:
        return config["variant_path"]
    else:
        return f"c7_graph_construction/graph_merge/{config['prefix']}.{{chr}}.variant.path"
        
def get_hap1_fa_input(wildcards):
    if "hap1_adaptor_masked_fa" in config:
        return config["hap1_fa"]
    else:
        return "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap1.fasta"

def get_hap2_fa_input(wildcards):
    if "hap2_adaptor_masked_fa" in config:
        return config["hap2_fa"]
    else:
        return "c6_draft_assembly/sample_assembly/{sample}/assembly/{sample}.hap2.fasta"
