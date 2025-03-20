# Pangenome construction and simplification
This pipeline performs pangenome construction and simplification.

<img align="middle" width="740" src="pangenome_construction.jpg"/>

## split_minigraph
### Description
-  Construct the SV pangenome graph with minigraph.
-  Align assemblies to pangenome graph with minigraph and filter the low-quality aligments.
-  Divide the pangenome graph into subgraphs.
-  Split alignments and assembly sequences into subgraphs.
### Requirement
-  minigraph
-  seqkit
-  VG
-  odgi
-  Cactus
-  cactus-gfa-tools
-  Samtools
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/simplify_ml_pangenome.yaml
```

## graph_construction
### Description
-  Induce base-level pangenome graph with seqwish.
-  Smooth the pangenome graph with smoothxg.
-  Normalize the pangenome graph with GFAffix.
-  Clip the dna-brnn masked or minigraph unaligned regions.
-  Reduce the graph complexity by untangling multiple-copy nodes among haplotypes.
### Requirement
-  seqiwsh
-  smoothxg
-  VG
-  GFAffix  
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/graph_construction.yaml
```

## simplify_ml_pangenome
### Description
-  Select the training and testing nodes/edges based on training and testing samples.
-  Generate pangenome graph without training and testing samples.
-  Extracting the node/edge features, including allele count-related, node sequence-related and graph motif-related features. 
-  Training the Multi-Layer Perceptron (MLP) model to distinguish error and variant nodes/edges.
-  Filter the error node/edges and compact the pangenome graph.
-  Inject SNV haplotypes into pangenome as paths.
### Requirement
-  PyTorch
-  VG
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/simplify_ml_pangenome.yaml
```

## merge_pangenome
### Description
-  Merge the subgraph pangenome into the final pangenome graph.
### Requirement
-  VG
### Usage
```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/merge_pangenome.yaml
```
