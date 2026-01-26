# SNV detection

This pipeline performs population-level SNV detection by leveraging PacBio long reads and Illumina short reads.

<img align="middle" width="800" src="snv_detection.jpg"/>

## call_sr_snv

### Description

- Detect SNVs using Illumina reads with GATK.

### Requirement

- [BWA](https://github.com/lh3/bwa)
- [Samtools](https://github.com/samtools/samtools)
- [GATK](https://github.com/broadinstitute/gatk)
- [Bcftools](https://github.com/samtools/bcftools)

### Configuration(call_sr_snv.yaml)

The configuration file should contain:

**`samples`**: Specify the path to a text file listing all samples.
The file must be **space-delimited**, with:

- **Column 1:** Sample name
- **Column 2:** Sample sex

**`sr_fastqs`**: Specify the paths to the **paired-end short-read FASTQ files**.
Use `{sample}` as a wildcard, which will be replaced automatically by the sample names defined in `samples`.

**`reference`**:

- `CHM13`: Provide the path to the **T2T CHM13 human reference genome** (FASTA format, indexed by `samtools faidx` and `bwa index`).

**`GATK_Resource`**: Specify paths to the required **GATK reference resource datasets**, all lifted over to CHM13 coordinates:

- `hapmap` – HapMap 3.3
- `omni` – 1000G Omni 2.5
- `1000G` – 1000 Genomes phase 1 high-confidence SNPs
- `known_indel` – Known indels in hg38 assembly.
- `mills` – Mills + 1000G gold-standard indels
- `axiomPoly` – Axiom Exome Plus polymorphism dataset

**`prefix`**: Prefix used for naming output files.

### Usage

```bash
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/call_sr_snv.yaml --workflow-profile ./profile/config_slurm/
```

### Output
The output file should contain:
- **`c1_call_sr_snv/merged_vcf/{config['prefix']}.gatk.variant_recalibrated.filter.vcf.gz`**: Short-read SNV callset

## call_lr_snv

### Description

- Detect candidate SNVs using PacBio HiFi reads with DeepVariant.
- Jointly call SNVs at the population level with GLnexus.
- Genotype SNVs using PacBio ZMW reads with PacBio ZMW reads with WhatsHap.
- Refine the SNV genotypes from genotype likelihood with Beagle4.

### Requirement

- [minimap2](https://github.com/lh3/minimap2)
- [DeepVariant](https://github.com/google/deepvariant)
- [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- [Bcftools](https://github.com/samtools/bcftools)
- [WhatsHap](https://github.com/whatshap/whatshap)
- [Beagle4](https://faculty.washington.edu/browning/beagle/b4_1.html)

### Configuration(call_lr_snv.yaml)

The configuration file should contain:

**`samples`**: Specify the path to a text file listing all samples.
The file must be **space-delimited**, with:

- **Column 1:** Sample name
- **Column 2:** Sample sex

**`lr_hifi_fastqs`**: Specify the path(s) to **PacBio HiFi long-read FASTQ files** for each sample.
Use `{sample}` as a wildcard; it will be automatically replaced with sample names from the `samples` file.

**`lr_zmw_fastqs`**: Specify the path(s) to **PacBio ZMW FASTQ files**.
Use `{sample}` as a wildcard; it will be automatically replaced with sample names from the `samples` file.

**`reference`**:

- `CHM13`: Provide the path to the **T2T CHM13 human reference genome** (FASTA format, indexed by `samtools faidx`).

**`prefix`**: Prefix used for naming output files.

### Usage

```bash
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/call_lr_snv.yaml --workflow-profile ./profile/config_slurm/
```

### Output
The output file should contain:
- **`c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.beagle.vcf.gz`**: Long-read SNV callset


## merge_snv

### Description

- Merge the Long-read (LR) and Short-read (SR) SNV callsets and select the SNV genotypes based on Hardy–Weinberg equilibrium P-value.
- Filter the SNVs for each individual based on kmers from short reads and HiFi reads with Merfin.

### Requirement

- [Bcftools](https://github.com/samtools/bcftools)
- [Merfin](https://github.com/arangrhie/merfin)

### Configuration(merge_snv.yaml)

The configuration file should contain:

**`sr_fastqs`**: Specify the paths to the **paired-end short-read FASTQ files**.
Use `{sample}` as a wildcard; it will be automatically replaced with sample names from the `samples` file.

**`lr_hifi_fastqs`**: Specify the path(s) to **PacBio HiFi long-read FASTQ files** for each sample.
Use `{sample}` as a wildcard; it will be automatically replaced with sample names from the `samples` file.

**`reference`**:

- `CHM13`: Provide the path to the **T2T CHM13 human reference genome** (FASTA format, indexed by `samtools faidx`).

**`prefix`**: Prefix used for naming output files.

**`sr_vcf`** _(optional)_: Path to the short-read SNV VCF file.
If omitted, the workflow will use the default:
`c1_call_sr_snv/merged_vcf/{prefix}.gatk.variant_recalibrated.filter.analysis_set.biallelic.vcf.gz`

**`lr_vcf`** _(optional)_: Path to the long-read SNV VCF file.
If omitted, the workflow will use the default:
`c2_call_lr_snv/lr_beagle/concat/{prefix}.deepvariant.whatshap.filter.analysis_set.biallelic_snp.beagle.filter.vcf.gz`

### Usage

```bash
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/merge_snv.yaml --workflow-profile ./profile/config_slurm/
```

### Output
The output file should contain:
- **`c3_merge_snv/merged_vcf/{config['prefix']}.consensus.merfin.vcf.gz`**: SNV callset combining short-read and long-read callsets