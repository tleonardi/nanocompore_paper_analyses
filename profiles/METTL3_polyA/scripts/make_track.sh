#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
RESULTS="${BASEDIR}/profiles/METTL3_polyA/results"
BED="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_sig_sites_GMM_logit_pvalue_thr_0.05.bed"
BEDG="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_sig_sites_GMM_logit_pvalue_thr_0.05.bedgraph"

(head -n 1 $BED
tail -n +2 $BED | bedparse convertChr --assembly hg38 --target ucsc 
)> $RESULTS/out_sig_sites_GMM_logit_pvalue_thr_0.05.ucsc.bed

(head -n 1 $BEDG | sed 's/Sites/Sites - Signal/'
tail -n +2 $BEDG | bedparse convertChr --assembly hg38 --target ucsc | sort -k1,1 -k2,2n -k4,4nr | sort -k1,1 -k2,2n -u
)> $RESULTS/out_sig_sites_GMM_logit_pvalue_thr_0.05.ucsc.bedgraph
