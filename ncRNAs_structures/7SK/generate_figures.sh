#!/bin/bash

# RF00100.stockholm.txt was edited by hand to replace the gap at position 82 with a G, which is actually
# present in hg38

BASEDIR=$(git rev-parse --show-toplevel)
WD="${BASEDIR}/ncRNAs_structures/7SK/"
nanocomp_res_kd="${BASEDIR}/nanocompore_pipelines/METTL3_KD_ncRNAs/results/nanocompore/out_nanocompore_results.tsv"
#nanocomp_res_kd="${BASEDIR}/7sk/data/KD/nanocompore/out_nanocompore_results.tsv"
nanocomp_res_ivt="${BASEDIR}/nanocompore_pipelines/IVT_7SK/results/nanocompore/out_nanocompore_results.tsv"
#nanocomp_res_dkc="${BASEDIR}/nanocompore_pipelines/DKC1_ncRNAs/results/nanocompore/out_nanocompore_results.tsv"
nanocomp_res_dkc="${BASEDIR}/7sk/data/DKC1/nanocompore/out_nanocompore_results.tsv"
nanocomp_bed="${BASEDIR}/7sk/data/IVT/references/reference_transcriptome.bed"                                                                                                                                              
template="${WD}/7sk.sto" 
tx="ENST00000636484"
thr="0.05"
pval_col=8
real_start=52995620

r2r  --GSC-weighted-consensus ${WD}/RF00100.stockholm.txt $template 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1
cat <(head -n -1 $template) ${WD}/7sk_custom_annots.txt <(echo //) > ${WD}/7sk_consensus_structure.sto
cat <(head -n -1 $template) ${WD}/7sk_custom_annots.txt <(python ${WD}/../create_annotations.py $nanocomp_res_kd $nanocomp_bed $template $tx 0.01 8 $real_start) <(echo //) > ${WD}/7sk_mettl3.sto

cat <(head -n -1 $template) ${WD}/7sk_custom_annots.txt <(python ${WD}/../create_annotations.py $nanocomp_res_ivt $nanocomp_bed $template $tx 0.01 6 $real_start) <(echo //) > ${WD}/7sk_ivt.sto

cat <(head -n -1 $template) ${WD}/7sk_custom_annots.txt <(python ${WD}/../create_annotations.py $nanocomp_res_dkc $nanocomp_bed $template $tx 0.1 8 $real_start) <(echo //) > ${WD}/7sk_dkc1.sto

r2r ${WD}/7sk.meta ${WD}/7sk.pdf
