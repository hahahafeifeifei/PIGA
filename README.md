# Pangenome-Informed Genome Assembly (PIGA)

This is a workflow called Pangenome-Informed Genome Assembly (PIGA) for **population-scale diploid genome assembly**, which leverages pangenome graph to improve accuracy and resolve complex structural variations.



***TODO***

- successfully run the whole pipeline with **test files** (currently dry run is ok)
- finish rules of `draft assembly` (preparation of the personal variants set)
- How to generate the files in `some_files` directory?
- Refine the **name** of each rule,  add the **annotation** for each rule, maybe edit the **path** of files...

## Features
0_0

## Install

Please make sure the tools listed below have been installed and executable:
- A
- B
- C

## Run Test Example
Currently ony test example is ready to run.
```bash
snakemake -s Snakefile --cores 4 --configfile config.yaml
```

## Full Workflow Execution
1. Configure parameters in `config.yaml`:
     `samples`: sample list of the data.

     `training_samples`:

     `test_samples`:

     `sr_fastqs`: the path of short read sequencing data of each sample.

     `lr_HIFI_fastqs`: the path of hifi read sequencing data of each sample.

     `lr_ZMW_fastqs`: the path of zmw read sequencing data of each sample.

     `reference`: the path of human reference genome(`GRCh38` and `CHM13`)

     `GATK_Resource`: the path of reference panel from dbsnp,hapmap,omni,1000G and so on which would be useful in the GATK calling process.

     `prefix`: (default: `1kcp`) 

3. Run:

    ```bash
    snakemake --cores 64 --configfile config/config.yaml
    ```
## Input Requirements

| **Input data**     | Format    | Description                           |
| ------------------ | --------- | ------------------------------------- |
| Short Reads (SR)   | `FASTQ`   | Illumina/WGS reads (gzip required)    |
| Long Reads (LR)    | `FASTQ`   | PacBio HiFi/ZMW reads (gzip required) |
| Outgroup Pangenome | Minigraph | Pre-built pangenome graph             |

## Module Descriptions



| Rule File                     | Inputs                                                       | Outputs                                     | Description |
| ----------------------------- | ------------------------------------------------------------ | ------------------------------------------- | ----------- |
| `call_sr_snv`                 | population SR read                                           | population SR snv                           |             |
| `call_lr_snv`                 | population LR read                                           | population LR snv                           |             |
| `merge_snv`                   | population SR snv + population LR snv                        | population merge snv                        |             |
| `phase_snv`                   | population merge snv + population LR read                    | population phased snv                       |             |
| `generate_personal_reference` | individual SR + individual LR + outgroup pangenome +individual phased snv | personal reference + personal reference snv |             |
| *To be finished*              |                                                              |                                             |             |
| `split_minigraph`             | Population FASTA                                             | Split subgraph chunks                       |             |
| `construct_pangenome`         | Split subgraph chunks                                        | Split subgraph pangenome                    |             |
| `simplify_ml_pangenome`       | Split subgraph pangenome                                     | Simplified subgraph pangenome               |             |
| `inject_snv_pangenome`        | Simplified subgraph pangenome                                | injected subgraph pangenome                 |             |
| `merge_pangenome`             | all iniected subgraph pangenomes                             | merged pangenome                            |             |
| `infer_diploid_path`          | individual SR + individual LR + individual draft assembly + merged pangenome | individual final assembly                   |             |



## Citation



## License

**License**: MIT License

