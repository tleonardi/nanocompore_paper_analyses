#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
ANALYSIS="${BASEDIR}/profiles/METTL3_polyA/analysis"
RESULTS="${BASEDIR}/profiles/METTL3_polyA/results"
NANOCOMPORE_FULL="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv"

mkdir -p $BASEDIR

cut -f 2,3,6,9 $NANOCOMPORE_FULL | awk 'BEGIN{OFS=FS="\t"}NR>1 && $4<0.5{$4=-log($4)/log(10); print $1,$2,$3,$4}' | sort -k1,1 -k2,2 -k4,4gr | sort -k1,1 -k2,2n -u | sort -k4,4gr | awk '{print ">"$1"_"$2"_"$4"\n"$3}' > $ANALYSIS/motifs/full_sites_ext.fa
grep ">" $ANALYSIS/motifs/full_sites_ext.fa | tr -d '>' > $ANALYSIS/motifs/full_sites_ext.fa.ids.txt                                                                                                                                                                                          

sylamer -fasta $ANALYSIS/motifs/full_sites_ext.fa -universe $ANALYSIS/motifs/full_sites_ext.fa.ids.txt -k 5 -grow 100 --over  -o $ANALYSIS/motifs/sylamer_full.out

