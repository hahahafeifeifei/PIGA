#!/bin/bash
# Download the required test dataset

# Sequencing data of test samples
wget -c -i test_sample.url.txt
for sample in {test-sample1,test-sample2,test-sample3}
do
    mkdir ${sample}
    mv ${sample}.* ${sample}
done

# GATK resource
wget -c -i gatk_resource.url.txt
mkdir reference/gatk_bundle
mv *vcf* reference/gatk_bundle

# Reference genome
## GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna reference
samtools index reference/15_GRCh38_no_alt_analysis_set
## CHM13
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz
mv chm13v2.0.fa reference
samtools index reference/chm13v2.0.fa

# External pangenome (HGSVC3)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Graph_Genomes/1.0/2024_02_23_minigraph_cactus_hgsvc3_hprc/hgsvc3-hprc-2024-02-23-mc-chm13.gbz
mv hgsvc3-hprc-2024-02-23-mc-chm13.gbz reference