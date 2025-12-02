# Tutorial

## Full Workflow Execution

```bash
# slurm system
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/config.yaml --workflow-profile ./profile/config_slurm/
# local run
snakemake -s Snakefile --cores 64 --jobs 64 --configfile config/config.yaml --workflow-profile ./profile/config_local/
```

## Running step by step

- [SNV detection](01_snv_detection.md)
- [SNV phasing](02_snv_phasing.md)
- [Personalized reference generation](03_personal_reference.md)
- [Draft diploid assembly](04_draft_assembly.md)
- [Pangenome construction and simplification](05_pangenome_construction.md)
- [Final diploid assembly reconstruction](06_final_assembly.md)
