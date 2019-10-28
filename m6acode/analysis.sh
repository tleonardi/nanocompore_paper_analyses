#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)

conda activate $BASEDIR/m6acode

python generate_new_db_actin.py
python generate_new_db_7sk.py

python parse_sampcompdb.py out_actin $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/references/reference_transcriptome.fa ENST00000331789 [1533,650,1322]
python parse_sampcompdb.py out_7sk $BASEDIR/nanocompore_pipelines/METTL3_KD_ncRNAs/results/references/reference_transcriptome.fa ENST00000636484 all


