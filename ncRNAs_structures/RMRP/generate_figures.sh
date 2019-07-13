#!/bin/bash

# RF00100.stockholm.txt was edited by hand to replace the gap at position 82 with a G, which is actually
# present in hg38

BASEDIR=$(git rev-parse --show-toplevel)
WD="${BASEDIR}/ncRNAs_structures/RMRP/"
nanocomp_res_kd="${BASEDIR}/nanocompore_pipelines/METTL3_KD_ncRNAs/results/nanocompore/out_nanocompore_results.tsv"
nanocomp_res_ivt="${BASEDIR}/nanocompore_pipelines/IVT_7SK/results/nanocompore/out_nanocompore_results.tsv"
nanocomp_res_dkc="${BASEDIR}/nanocompore_pipelines/DKC1_ncRNAs/results/nanocompore/out_nanocompore_results.tsv"
#nanocomp_res_dkc="${BASEDIR}/7sk/data/DKC1/nanocompore/out_nanocompore_results.tsv"
nanocomp_bed="${BASEDIR}/7sk/data/IVT/references/reference_transcriptome.bed"

template="${WD}/rmrp.sto" 
tx="ENST00000602361"
thr="0.05"
pval_col=8
real_start=35657752
rfam="M29212.1/2-265"

r2r  --GSC-weighted-consensus ${WD}/RF00030.stockholm.txt $template 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1
cat <(head -n -1 $template) ${WD}/rmrp_custom_annots.txt <(echo //) > ${WD}/rmrp_consensus_structure.sto
#cat <(head -n -1 $template) ${WD}/rmrp_custom_annots.txt <(python ${WD}/../create_annotations.py $nanocomp_res_kd $nanocomp_bed $template $tx 0.01 8 $real_start $rfam) <(echo //) > ${WD}/rmrp_mettl3.sto

cat <(head -n -1 $template) ${WD}/rmrp_custom_annots.txt <(python ${WD}/../create_annotations.py $nanocomp_res_dkc $nanocomp_bed $template $tx 0.01 8 $real_start $rfam) <(echo //) > ${WD}/rmrp_dkc1.sto

r2r ${WD}/rmrp.meta ${WD}/rmrp.pdf
