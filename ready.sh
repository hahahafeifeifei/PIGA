#!/bin/sh
cp Snakefile.py Snakefile
cp rules/SR_var_calling.py rules/SR_var_calling.smk
cp rules/LR_var_calling.py rules/LR_var_calling.smk
cp rules/merge_snv.py rules/merge_snv.smk
cp rules/phase_snv.py rules/phase_snv.smk
cp rules/draft_assembly.py rules/draft_assembly.smk
cp rules/split_minigraph.py rules/split_minigraph.smk
cp rules/graph_construction.py rules/graph_construction.smk
cp rules/simplify_ml_pangenome.py rules/simplify_ml_pangenome.smk
cp rules/merge_pangenome.py rules/merge_pangenome.smk
cp rules/infer_diploid_path.py rules/infer_diploid_path.smk
cp rules/generate_personal_reference.py rules/generate_personal_reference.smk
