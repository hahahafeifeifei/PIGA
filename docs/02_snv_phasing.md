# SNV phasing pipeline
This pipeline performs **SNV phasing** by integrating information from long reads and population.

<img align="middle" width="800" src="snv_phasing.jpg"/>

## phase_snv
### Description
-  Generate long read-based SNV phasing blocks for each individual with WhatsHap.
-  Perform statistical phasing with Shapeit4 by intergrating the prior haplotype information from long read and external variant panel.
### Requirement
-  Bcftools
-  WhatsHap
-  Shapeit4
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/phase_snv.yaml
```
