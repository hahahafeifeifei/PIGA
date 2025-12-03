# Draft diploid genome assembly

This pipeline performs reference-guided draft diploid genome assembly.

<img align="middle" width="942" src="draft_assembly.jpg"/>

## draft_assembly

### Description

- Partition long reads into haplotypes with WhatsHap.
- Detect intervals with enough partitioned long-read coverage and extract their contained reads to perform assembly.
- Trim the subreads with MECAT2.
- Perform subinterval assembly with wtdbg2.
- Polish contigs using subreads with Racon and Arrow.
- Recalibrate the base quality of contigs and filter low-quality contigs.
- Perform interval assembly with phasebook.
- Detect structural errors in contigs with Sniffles2 and Jasmine.
- Polish contigs using short reads with Freebayes.
- Filter contig regions with switch errors by comparing assembly haplotype to SNV haplotype.
- Mask the adaptor sequences in contigs.

### Requirement

- [BWA](https://github.com/lh3/bwa)
- [GATK](https://github.com/broadinstitute/gatk)
- [minimap2](https://github.com/lh3/minimap2)
- [WhatsHap](https://github.com/whatshap/whatshap)
- [Samtools](https://github.com/samtools/samtools)
- [Bcftools](https://github.com/samtools/bcftools)
- [Seqkit](https://github.com/shenwei356/seqkit)
- [MECAT2](https://github.com/xiaochuanle/MECAT2)
- [wtdbg2](https://github.com/ruanjue/wtdbg2)
- [Racon](https://github.com/isovic/racon)
- [Arrow](https://github.com/PacificBiosciences/gcpp#)
- [cutadapt](https://github.com/marcelm/cutadapt/)
- [phasebook](https://github.com/phasebook/phasebook)
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [Jasmine](https://github.com/mkirsche/Jasmine)
- [Freebayes](https://github.com/freebayes/freebayes)
- [Bedtools](https://github.com/arq5x/bedtools2)

### Configuration (`draft_assembly.yaml`)

The configuration file should contain:

**`samples`**: Specify the path to a text file listing all samples.
The file must be **space-delimited**, with:

- **Column 1:** Sample name
- **Column 2:** Sample sex

**`sr_fastqs`**: Specify the paths to the **paired-end short-read FASTQ files**.
Use `{sample}` as a wildcard; it will be automatically replaced with sample names from the `samples` file.

**`lr_hifi_fastqs`**: Specify the path(s) to **PacBio HiFi long-read FASTQ files**.
Use `{sample}` as a wildcard.

**`lr_zmw_fastqs`**: Specify the path(s) to **PacBio ZMW FASTQ files**.
Use `{sample}` as a wildcard.

**`lr_subreads_bam`**: Specify the path(s) to **PacBio subreads BAM files**.
Use `{sample}` as a wildcard.

**`reference`**:

- `CHM13`: Provide the path to the **T2T CHM13 human reference genome** (FASTA format, indexed by `samtools faidx`).

**`par_region`**: Path to the **pseudoautosomal region (PAR) BED file** for CHM13.
This is used to correctly handle sex-chromosome reconstruction.

**`prefix`**: Prefix used for naming output files.

**`consensus_fasta`** _(optional)_: Path to the sample-specific consensus FASTA file.
If omitted, the workflow will use the default:
`c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.fasta`

**`chain`** _(optional)_: Path to the chain file describing coordinate mapping between CHM13 and the consensus assembly.
If omitted, the workflow will use the default:
`c5_personal_ref/consensus_fasta/{sample}/CHM13.af_pangenome.{sample}_polish.chain`

**`sample_phased_vcf`** _(optional)_: Path to the sample's phased VCF.
If omitted, the workflow will use the default:
`c4_phase_snv/shapeit4/samples/{sample}/{sample}.shapeit4.vcf.gz`

**`sample_WGS_meryl`** _(optional)_: Path to the sample’s WGS meryl directory.
If omitted, the workflow will use the default:
`c3_merge_snv/meryl/{sample}/{sample}-WGS.meryl`

### Usage

```bash
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/draft_assembly.yaml --workflow-profile ./profile/config_slurm/
```
