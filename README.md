# Pangenome-Informed Genome Assembly (PIGA)

## Introduction

This is a workflow called Pangenome-Informed Genome Assembly (PIGA) for **population-scale diploid genome assembly**, which employs a pangenome graph as a unified framework to integrate sequence information across individuals and perform joint diploid genome assembly.

Compared to the current assembly methods, the PIGA workflow fully utilizes multiple sources of information (**long read**, **short read**, **internal population**, **external assembly panel**) and is well adapted to **low-coverage** and **modest-coverage** situations.

<img align="middle" width="1000" src="PIGA.jpg"/>

## Installation & Setup

#### Step 1: Clone the Repository

First, clone the PIGA repository from GitHub:

```bash
git clone https://github.com/JianYang-Lab/PIGA.git
cd PIGA
```

#### Step 2: Create and Activate Conda Environment

PIGA uses Snakemake and requires several software dependencies. The easiest way to install them is by creating a dedicated conda environment from the provided YAML file.

```bash
# Create the conda environment named 'piga'
conda env create -f environment.yaml -n piga

# Activate the environment to use PIGA
conda activate piga
```

#### Step 3: Build and Download Additional Tools

Some tools are not available via conda and need to be built from source or downloaded separately.

```bash
# Build and download additional tools
bash build.sh
```

## Test Dataset and Configuration

#### Download the test dataset

```bash
cd test_data
bash download.sh
```

#### Configure the workflow(config/config.yaml)

The configuration file should contain:

**`samples`**: Specify the path to a text file listing all samples.
The file must be **space-delimited**, with:

- **Column 1:** Sample name
- **Column 2:** Sample sex

**`sr_fastqs`**: Specify the paths to the **paired-end short-read FASTQ files**.
Use `{sample}` as a wildcard.

**`lr_hifi_fastqs`**: Specify the path(s) to **PacBio HiFi long-read FASTQ files**.
Use `{sample}` as a wildcard.

**`lr_zmw_fastqs`**: Specify the path(s) to **PacBio ZMW FASTQ files**.(Selected representative read for each ZMW)
Use `{sample}` as a wildcard.

**`lr_subreads_bam`**: Specify the path(s) to **PacBio subreads BAM files**.
Use `{sample}` as a wildcard.

**`reference`**:

- `CHM13`: Path to the **T2T CHM13 reference genome** (FASTA format, indexed by `samtools faidx` and `bwa index`).
- `GRCh38`: Path to the **GRCh38 reference genome**in FASTA format.

**`GATK_Resource`**: Paths to **GATK reference resource datasets**:

- `hapmap` – HapMap 3.3

- `omni` – 1000G Omni 2.5

- `1000G` – 1000 Genomes phase 1 high-confidence SNPs

- `known_indel` – Known indels in hg38.

- `mills` – Mills + 1000G gold-standard indels

- `axiomPoly` – Axiom Exome Plus polymorphism dataset

**`external_pangenome`**: Path to the **external pangenome graph** (GBZ format).

**`par_region`**: Path to the **pseudoautosomal region (PAR) BED file** for CHM13.

**`internal_assembly_list`**: Path to a file listing **internal genome assemblies**.

**`external_assembly_list`**: Path to a file listing **external genome assemblies**.

**`train_sample_list`**: Path to a file listing samples used for model training.

**`prefix`**: Prefix used for naming output files.
You are now ready to run the workflow.

**TODO**: tell users how to set up the configuration file including config/config.yaml and also other step-by-step config files.

## Running the Workflow

After setting up the configuration file (`config/config.yaml`) with your sample data, you can run the entire PIGA workflow. There are two main ways to execute it:

#### Option 1: Local Execution

This method is suitable for running PIGA on a single, powerful machine. It will use the specified number of cores on the local machine.

```bash
# Run PIGA locally using up to 64 cores
snakemake -s Snakefile --cores 64 --configfile config/config.yaml --workflow-profile ./profile/config_local.yaml
```

#### Option 2: Cluster Execution (with Profile)

This method is designed for high-performance computing (HPC) environments and uses a job submission system (e.g., SLURM) to distribute the workload across a cluster.

```bash
# Run PIGA using a workflow profile to submit jobs to a cluster
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/config.yaml --workflow-profile ./profile/config_slurm.yaml
```

**Note:** By default, we provide a profile which is configured to use the **SLURM** job scheduler. You can customize the cluster settings (e.g., switch to a different scheduler or change resource allocation) by editing the configuration file.

## Documentation

PIGA consists of six modules, each containing several commands that can be executed independently. A detailed tutorial is provided for each module.

#### [1. SNV Detection](docs/01_snv_detection.md)

- `call_sr_snv`: detect SNVs using short reads.
- `call_lr_snv`: detect SNVs using long reads.
- `merge_snv`: merge the short-read SNV callset and long-read SNV callset.

#### [2. SNV Haplotype Phasing](docs/02_snv_phasing.md)

- `phase_snv`: perform SNV haplotype phasing leveraging long-read and population information.

#### [3. Personalized Reference Generation](docs/03_personal_reference.md)

- `generate_personal_reference`: generate personalized reference by modifying the reference genome with homozygous variants genotyped from the external pangenome.

#### [4. Draft Diploid Genome Assembly](docs/04_draft_assembly.md)

- `draft_assembly`: partition long reads into haplotypes and produce draft diploid assemblies.

#### [5. Pangenome Construction and Simplification](docs/05_pangenome_construction.md)

- `construct_pangenome`: construct and refine the base-level pangenome.
- `simplify_pangenome`: simplify the pangenome.
- `merge_pangenome`: merge pangenome subgraphs into the final pangenome.

#### [6. Final Diploid Assembly Reconstruction](docs/06_final_assembly.md)

- `infer_diploid_path`: reconstruct the final diploid assembly by inferring the diploid paths.

## License

**License**: MIT License
