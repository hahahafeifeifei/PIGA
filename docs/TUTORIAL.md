# TUTORIAL

## Full Workflow Execution

```shell
snakemake -s Snakefile --cores 64 --configfile config.yaml --configfile tools.yaml
```



## Running only one step



**(1) call_sr_snv**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/SR_var_calling.yaml
```



**(2) call_lr_snv**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/LR_var_calling.yaml
```



**(3) merge_snv**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/merge_snv.yaml
```



**(4) phase_snv**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/phase_snv.yaml
```



**(5) generate_personal_reference**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/generate_personal_reference.yaml
```



**(6) draft_assembly**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/draft_assembly.yaml
```



**(7) split_minigraph**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/split_minigraph.yaml
```



**(8) graph_construction**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/graph_construction.yaml
```



**(9) simplify_ml_pangenome**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/simplify_ml_pangenome.yaml
```



**(10) merge_pangenome**

```shell
snakemake -s Snakefile --cores 64 --configfile config/tools.yaml --configfile config/merge_pangenome.yaml
```



**(11) infer_diploid_path**

```shell
snakemake -s Snakefile --cores 64 --configfile config/config.yaml --configfile config/tools.yaml --configfile config/infer_diploid_path.yaml
```



