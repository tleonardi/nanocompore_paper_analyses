#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)


DATASETS="METTL3_KD_polyA TRM5_KO_polyA"
for dataset in $DATASETS; do
	mkdir -p $BASEDIR/generic_seq_metrics/results/$dataset
	(
	for f in $BASEDIR/nanocompore_pipelines/$dataset/results/*/minimap.filt.sort.bam; do
	echo -ne "$(basename $(dirname $f))\t"
	samtools view $f | cut -f1 | sort -u | wc -l
	done
	)> $BASEDIR/generic_seq_metrics/results/$dataset/n_mapped_reads.txt
	
	(
	for f in $BASEDIR/nanocompore_pipelines/$dataset/results/*/minimap.filt.sort.bam; do
	SAMPLE=$(basename $(dirname $f))
	samtools view $f | cut -f3 | sort | uniq -c | awk -v samp=$SAMPLE 'BEGIN{OFS="\t"}{print samp, $1,$2}' 
 	done
	)> $BASEDIR/generic_seq_metrics/results/$dataset/tx_counts.txt

done

dataset="DKC1_ncRNAs"
mkdir -p $BASEDIR/generic_seq_metrics/results/$dataset
(
for f in $BASEDIR/nanocompore_pipelines/$dataset/results/*/minimap.filt.sort.bam; do
SAMPLE=$(basename $(dirname $f))
samtools view $f | cut -f3 | sort | uniq -c | awk -v samp=$SAMPLE 'BEGIN{OFS="\t"}{print samp, $1,$2}' 
done
)> $BASEDIR/generic_seq_metrics/results/$dataset/tx_counts.txt
