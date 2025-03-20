# SNV phasing
This pipeline performs SNV phasing by integrating information from long reads and population.

<img align="middle" width="645" src="snv_phasing.jpg"/>

## phase_snv
### Description
-  Generate long read-based SNV phasing blocks for each individual with WhatsHap.
-  Perform statistical phasing by intergrating the prior haplotype information from long read and external haplotype panel with Shapeit4.
### Requirement
-  Bcftools
-  WhatsHap
-  Shapeit4
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/phase_snv.yaml
```
