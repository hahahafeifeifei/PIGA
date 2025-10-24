def get_sr_input_fastqs(wildcards):
    return config["sr_fastqs"][wildcards.sample]

def get_zmw_input_fastqs(wildcards):
    return config["lr_zmw_fastqs"][wildcards.sample]

def get_hifi_input_fastqs(wildcards):
    return config["lr_hifi_fastqs"][wildcards.sample]

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
        return f"c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.analysis_set.biallelic.vcf.gz"

def get_lr_vcf_input(wildcards):
    if "lr_vcf" in config:
        return config["lr_vcf"]
    else:
        return f"c2_call_lr_snv/lr_beagle/concat/CKCG.deepvariant.whatshap.filter.analysis_set.biallelic_snp.beagle.filter.vcf.gz"

def get_merge_merfin_filter_vcf_input(wildcards):
    if "merge_merfin_filter_vcf" in config:
        return config["merge_merfin_filter_vcf"]
    else:
        return "c3_merge_snv/merfin/merge/CKCG.CHM13.consensus.phase1.call_set.hwe_missing_filter.vcf.gz"

def get_sr_scaffold_vcf_input(wildcards):
    if "sr_scaffold_vcf" in config:
        return config["sr_scaffold_vcf"]
    else:
        return "c3_merge_snv/callset/srs_lrs_compare/CKCG.analysis_set.call_set.consensus.srs_scaffold.vcf.gz"

def get_zmw_bam_input(wildcards):
    if "zmw_bam" in config:
        return config["zmw_bam"]
    else:
        return "c2_call_lr_snv/lr_mapping/{sample}/{sample}.zmw.pbmm2.bam"

def get_consensus_fasta_input(wildcards):
    if "consensus_fasta" in config:
        return config["consensus_fasta"]
    else:
        return "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta"

def get_chain_input(wildcards):
    if "chain" in config:
        return config["chain"]
    else:
        return "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.chain"

def get_sample_vcf_input(wildcards):
    if "sample_vcf" in config:
        return config["sample_vcf"]
    else:
        return "c4_phase_snv/shapeit4/samples/{sample}/{sample}.shapeit4.vcf.gz"

def get_sample_meryl_input(wildcards):
    if "sample_meryl" in config:
        return config["sample_meryl"]
    else:
        return "c3_merge_snv/meryl/{sample}/{sample}.meryl/merylIndex"

def concat_final_phase_vcf_sex_specific_chrlist(wildcards):
    sex = config['sex'][wildcards.sample]
    chr_list = [f'chr{i}' for i in range(1, 23)]
    if sex == 'female':
        return chr_list + ['chrX']
    return chr_list


def get_hapl_input(wildcards):
    if "hapl" in config:
        return config["hapl"]
    else:
        return "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.hapl"

def get_gbz_input(wildcards):
    if "gbz" in config:
        return config["gbz"]
    else:
        return "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.merge.assembly.gbz"

def get_gfa_input(wildcards):
    if "gfa" in config:
        return config["gfa"]
    else:
        return "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.assembly.gfa"

def get_variant_path_input(wildcards):
    if "variant_path" in config:
        return config["variant_path"]
    else:
        return "c7_graph_construction/graph_merge/t2t.grch38.58hifi.1064zmw.{chr}.variant.path"
        
def get_hap1_adaptor_masked_fa_input(wildcards):
    if "hap1_adaptor_masked_fa" in config:
        return config["hap1_adaptor_masked_fa"]
    else:
        return "c6_draft_assembly/{sample}/assembly/{sample}.hap1.adaptor_masked.fasta"

def get_hap2_adaptor_masked_fa_input(wildcards):
    if "hap2_adaptor_masked_fa" in config:
        return config["hap2_adaptor_masked_fa"]
    else:
        return "c6_draft_assembly/{sample}/assembly/{sample}.hap2.adaptor_masked.fasta"
