# Personalized reference generation
This pipeline generates personalized reference and SNV haplotype for each individual.

<img align="middle" width="800" src="personal_reference.jpg"/>

## generate_personal_reference
### Description
-  Genotype variants (SNVs, indels and SVs) from the external pangenome with VG and GraphAligner.
-  Generate personalized reference by modifying original reference using homozygous variants with Bcftools.
-  Lift SNV haplotypes from original reference to personalized reference.
-  Call SNVs for the lifted-failed regions with DeepVariant.
-  Phased newly called SNV into SNV haplotypes with WhatsHap.
### Requirement
-  VG
-  GraphAligner
-  Bcftool
-  GATK
-  BWA
-  DeepVariant
-  WhasHap
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/generate_personal_reference.yaml
```
