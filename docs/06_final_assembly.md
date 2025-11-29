# Final diploid assembly reconstruction
This pipeline reconstructs the final diploid assembly by inferring the diploid path from pangenome.

<img align="middle" width="667" src="final_assembly.jpg"/>

## infer_diploid_path
### Description
-  Sample the candidate haplotypes from pangenome based on kmers with VG.
-  Select the haplotype path with top score as haploid reference.
-  Normalize the personalized pangenome with VG and GFAffix.
-  Infer the diploid path by performing integrating genotyping and phasing using information from SNV path, draft assembly path, long read and short read.
-  Patch the diploid assembly by the draft assembly at structure-discordant regions.
-  Clip the diploid assembly at kmer-erroneous regions.
### Requirement
-  VG
-  GFAffix
-  vcfbub
-  Bcftools
-  Pangenie
-  MarginPhase
-  Hiphase
-  minimap2
-  Secphase
-  cuteSV
-  paftools
-  Merqury
-  Bedtools
-  Samtools
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/infer_diploid_path.yaml
```
