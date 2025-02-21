# rule all:
#     input:
#         expand("c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta", sample=config['samples'])

rule lr_GraphAligner_mapping:
    input:
        pbcc_fastq = "/storage/yangjianLab/wangyifei/project/01.CKCG/03.CCS+CLR/result/pbccs/{sample}-CLR/{sample}-CLR.pbccs_all.merge.fastq.gz"
    output:
        pbcc_gam = "c5_personal_ref/pbcc_gam/{sample}/{sample}-CLR.graphaligner.gam"
    threads: 16
    resources: 
        mem_mb = 200000
    params:
        pangenome_name = "/storage/yangjianLab/wangyifei/project/02.HanGenome/Pangenome/01.Cactus/v4/06.graphmap-join-clip.min12/t2t.grch38.62chineses.clip.min12"
    shell:
        """
        /storage/yangjianLab/duanzhongqu/miniconda3/bin/GraphAligner \
        -t {threads} \
        -g {params.pangenome_name}.vg \
        -f {input.pbcc_fastq} \
        -a {output.pbcc_gam} \
        -x vg
        """
rule sr_giraffe_mapping:
    input:
        ngs_fastq_1 = "/storage/yangjianLab/sharedata/CKCG/preprocess.data/WGS/fastp/{sample}-WGS/{sample}-WGS.fastp.R1.fastq.gz",
        ngs_fastq_2 = "/storage/yangjianLab/sharedata/CKCG/preprocess.data/WGS/fastp/{sample}-WGS/{sample}-WGS.fastp.R2.fastq.gz"
    output:
        ngs_gam = "c5_personal_ref/ngs_gam/{sample}/{sample}-WGS.giraffe.1.gam"
    threads: 16
    resources: 
        mem_mb = 200000
    params:
        pangenome_name = "/storage/yangjianLab/wangyifei/project/02.HanGenome/Pangenome/01.Cactus/v4/06.graphmap-join-clip.min12/t2t.grch38.62chineses.clip.min12"
    shell:
        """
        vg giraffe \
        -t {threads} \
        -H {params.pangenome_name}.gbwt \
        -g {params.pangenome_name}.gg \
        -m {params.pangenome_name}.min \
        -d {params.pangenome_name}.dist \
        -f {input.ngs_fastq_1} \
        -f {input.ngs_fastq_2} \
        > {output.ngs_gam}
        """
        
rule sr_giraffe_R1_unpair_mapping:
    input:
        ngs_unpair_fastq_1 = "/storage/yangjianLab/sharedata/CKCG/preprocess.data/WGS/fastp/{sample}-WGS/{sample}-WGS.fastp.unpaired.R1.fastq.gz"
    output:
        ngs_gam = "c5_personal_ref/ngs_gam/{sample}/{sample}-WGS.giraffe.2.gam"
    threads: 16
    resources: 
        mem_mb = 200000
    params:
        pangenome_name = "/storage/yangjianLab/wangyifei/project/02.HanGenome/Pangenome/01.Cactus/v4/06.graphmap-join-clip.min12/t2t.grch38.62chineses.clip.min12"
    shell:
        """
        vg giraffe \
        -t {threads} \
        -H {params.pangenome_name}.gbwt \
        -g {params.pangenome_name}.gg \
        -m {params.pangenome_name}.min \
        -d {params.pangenome_name}.dist \
        -f {input.ngs_unpair_fastq_1} \
        > {output.ngs_gam}
        """
        
rule sr_giraffe_R2_unpair_mapping:
    input:
        ngs_unpair_fastq_2 = "/storage/yangjianLab/sharedata/CKCG/preprocess.data/WGS/fastp/{sample}-WGS/{sample}-WGS.fastp.unpaired.R2.fastq.gz"
    output:
        ngs_gam = "c5_personal_ref/ngs_gam/{sample}/{sample}-WGS.giraffe.3.gam"
    threads: 16
    resources: 
        mem_mb = 200000
    params:
        pangenome_name = "/storage/yangjianLab/wangyifei/project/02.HanGenome/Pangenome/01.Cactus/v4/06.graphmap-join-clip.min12/t2t.grch38.62chineses.clip.min12"
    shell:
        """
        vg giraffe \
        -t {threads} \
        -H {params.pangenome_name}.gbwt \
        -g {params.pangenome_name}.gg \
        -m {params.pangenome_name}.min \
        -d {params.pangenome_name}.dist \
        -f {input.ngs_unpair_fastq_2} \
        > {output.ngs_gam}
        """
        
rule gam_merge:
    input:
        pbcc_gam = "c5_personal_ref/pbcc_gam/{sample}/{sample}-CLR.graphaligner.gam",
        ngs_gam_1 = "c5_personal_ref/ngs_gam/{sample}/{sample}-WGS.giraffe.1.gam",
        ngs_gam_2 = "c5_personal_ref/ngs_gam/{sample}/{sample}-WGS.giraffe.2.gam",
        ngs_gam_3 = "c5_personal_ref/ngs_gam/{sample}/{sample}-WGS.giraffe.3.gam"
    output:
        ngs_merged_gam = "c5_personal_ref/merged_gam/{sample}/{sample}.merge.gam"
    shell:
        """
        cat {input.pbcc_gam} {input.ngs_gam_1} {input.ngs_gam_2} {input.ngs_gam_3} > {output.ngs_merged_gam}
        
        rm {input.pbcc_gam}
        rm {input.ngs_gam_1}
        rm {input.ngs_gam_2}
        rm {input.ngs_gam_3}
        """
rule gam_call_variants:
    input:
        gam = "c5_personal_ref/merged_gam/{sample}/{sample}.merge.gam"
    output:
        vcf = "c5_personal_ref/vg_call/{sample}/{sample}.af_pangenome.merge.vcf"
    threads: 16
    resources:
        mem_mb = 200000
    params:
        pangenome_name = "/storage/yangjianLab/wangyifei/project/02.HanGenome/Pangenome/01.Cactus/v4/06.graphmap-join-clip.min12/t2t.grch38.62chineses.clip.min12"
    shell:
        """
        vg pack \
        -t {threads} \
        -x {params.pangenome_name}.vg \
        -g {input.gam} \
        -o c5_personal_ref/vg_call/{wildcards.sample}/{wildcards.sample}.merge.pack
        
        vg call \
        -t {threads} \
        -m 3,6 \
        {params.pangenome_name}.vg \
        -r {params.pangenome_name}.snarls \
        -g {params.pangenome_name}.gbwt \
        -k c5_personal_ref/vg_call/{wildcards.sample}/{wildcards.sample}.merge.pack \
        -s {wildcards.sample} \
        > {output.vcf}
        
        rm c5_personal_ref/vg_call/{wildcards.sample}/{wildcards.sample}.merge.pack
        """
        
rule gam_call_variants_filter:
    input:
        vcf = "c5_personal_ref/vg_call/{sample}/{sample}.af_pangenome.merge.vcf"
    output:
        filtered_vcf = "c5_personal_ref/vg_call/{sample}/{sample}.af_pangenome.merge.filter.vcf.gz"
    threads: 2
    shell:
        """
        dp_mean=$(bcftools query -f "%DP\n" {input.vcf} | awk '{{if($1>=100)print 100;else print$0}}' | csvtk -H -t summary -f 1:mean)                            
        dp_sd=$(bcftools query -f "%DP\n" {input.vcf} | awk '{{if($1>=100)print 100;else print$0}}' | csvtk -H -t summary -f 1:stdev)
        bcftools view -f PASS -i "GT='AA' && FORMAT/DP<=$(awk -v dp_sd=${{dp_sd}} -v dp_mean=${{dp_mean}} 'BEGIN{{print dp_mean+3*dp_sd}}' )" {input.vcf} | sed 's/T2Tv2.chr/chr/g' > c5_personal_ref/vg_call/{wildcards.sample}/{wildcards.sample}.af_pangenome.merge.filter.vcf
        
        bgzip -f -@ {threads} c5_personal_ref/vg_call/{wildcards.sample}/{wildcards.sample}.af_pangenome.merge.filter.vcf
        tabix -f {output.filtered_vcf}
        """

rule filtered_variants_ref_consensus:
    input:
        vcf = "c5_personal_ref/vg_call/{sample}/{sample}.af_pangenome.merge.filter.vcf.gz",
        ref = config['reference']
    output:
        consensus_fasta = "c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta"
    shell:
        """
        bcftools consensus -f {input.ref} -H 1 {input.vcf} > {output.consensus_fasta}
        """
        