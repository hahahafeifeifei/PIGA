# Pangenome-Informed Genome Assembly (PIGA)

## Introduction

This is a workflow called Pangenome-Informed Genome Assembly (PIGA) for **population-scale diploid genome assembly**, which employs a pangenome graph as a unified framework to integrate sequence information across individuals and perform joint diploid genome assembly.

Compared to the current assembly methods, the PIGA workflow fully utilizes multiple sources of information (long read, short read, internal population, external assembly panel) and is well adapted to low-coverage and modest-coverage situations.

<img align="middle" width="1000" src="PIGA.jpg"/>

## Installation & Setup

**Step 1: Clone the Repository**

First, clone the PIGA repository from GitHub:

```bash
git clone https://github.com/JianYang-Lab/PIGA.git
cd PIGA
```

**Step 2: Create and Activate Conda Environment**

PIGA uses Snakemake and requires several software dependencies. The easiest way to install them is by creating a dedicated conda environment from the provided YAML file.

```bash
# Create the conda environment named 'piga'
conda env create -f piga.yaml

# Activate the environment to use PIGA
conda activate piga
```

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

To be noted, the PIGA consists of several modules and alternatively you can run whole PIGA workflow step-by-step or just run some of the modules.

The following modules are available in PIGA:

- `call_sr_snv`: detect SNVs using short reads.
- `call_lr_snv`: detect SNVs using long reads.
- `merge_snv`: merge the short-read SNV callset and long-read SNV callset.
- `phase_snv`: perform SNV haplotype phasing leveraging long-read and population information.
- `generate_personal_reference`: generate personalized reference by modifying the reference genome with homozygous variants genotyped from the external pangenome.
- `draft_assembly`: partition long reads into haplotypes and produce draft diploid assemblies.
- `construct_pangenome`: construct and refine the base-level pangenome.
- `simplify_ml_pangenome`: simplify the pangenome.
- `merge_pangenome`: merge pangenome subgraphs into the final pangenome.
- `infer_diploid_path`: reconstruct the final diploid assembly by inferring the diploid paths.

The step-by-step documentation is available at [tutorial](docs/tutorial.md). If you encounter problems, please open a [github issue](https://github.com/JianYang-Lab/1KCP-analysis/issues).

## License

**License**: MIT License
