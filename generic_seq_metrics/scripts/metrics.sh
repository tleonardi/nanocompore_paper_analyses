#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)


mkdir -p $BASEDIR/generic_seq_metrics/results/METTL3_KD_polyA
(
echo -ne "WT_1\t"
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/WT_1/minimap.filt.sort.bam | cut -f1 | sort -u | wc -l
echo -ne "WT_2\t"
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/WT_2/minimap.filt.sort.bam | cut -f1 | sort -u | wc -l
echo -ne "KD_1\t"
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/KD_1/minimap.filt.sort.bam | cut -f1 | sort -u | wc -l
echo -ne "KD_2\t"
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/KD_2/minimap.filt.sort.bam | cut -f1 | sort -u | wc -l
)> $BASEDIR/generic_seq_metrics/results/METTL3_KD_polyA/n_mapped_reads.txt


(
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/WT_1/minimap.filt.sort.bam | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print "WT_1", $1,$2}' 
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/WT_2/minimap.filt.sort.bam | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print "WT_2", $1,$2}' 
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/KD_1/minimap.filt.sort.bam | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print "KD_1", $1,$2}' 
samtools view $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/KD_2/minimap.filt.sort.bam | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print "KD_2", $1,$2}' 
)> $BASEDIR/generic_seq_metrics/results/METTL3_KD_polyA/tx_counts.txt
