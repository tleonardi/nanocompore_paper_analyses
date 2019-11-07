#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
ANALYSIS="${BASEDIR}/profiles/METTL3_polyA/analysis/miCLIP"
data="$BASEDIR/profiles/METTL3_polyA/data/miCLIP_BAMs"
RESULTS="${BASEDIR}/profiles/METTL3_polyA/results"
NANOCOMPORE_FULL="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv"
REF_TRANSCRIPTOME="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/references/reference_transcriptome.bed"
GENOME="/home/tleonardi/nobackup/Homo_sapiens.GRCh38.dna.toplevel.fa"
TRANSCRIPTOME="/home/tleonardi/nobackup/Homo_sapiens.GRCh38.94.chr.bed"
conda activate $BASEDIR/conda_environments/star

mkdir -p $ANALYSIS/fastq
for bam in $data/*.bam; do
	name=$(basename $bam -jufastqgz_unmapped.bam.deduplicated.bam)
	samtools fastq -f 256 $bam > $analysis/fastq/${name}.fastq
done

(for bam in $data/*.bam; do
	name=$(basename $bam -jufastqgz_unmapped.bam.deduplicated.bam)
	echo -ne "$name\t"
	samtools view -F 256 $bam | wc -l 
done)>$ANALYSIS/tot_depth.txt

bedparse filter -a $NANOCOMPORE_FULL -c 4 $REF_TRANSCRIPTOME > $ANALYSIS/transcriptome.bed
bedtools getfasta -name -split -s -fi $GENOME -bed $ANALYSIS/transcriptome.bed -fo $ANALYSIS/transcriptome.fa
bedtools getfasta -name -split -s -fi $GENOME -bed $TRANSCRIPTOME -fo $ANALYSIS/whole_transcriptome.fa 


mkdir -p $ANALYSIS/bowties_indexes/
bowtie2-build $ANALYSIS/transcriptome.fa $ANALYSIS/bowties_indexes/transcriptome
mkdir -p $ANALYSIS/bowties_indexes_whole_transcriptome
bowtie2-build $ANALYSIS/whole_transcriptome.fa $ANALYSIS/bowties_indexes_whole_transcriptome/transcriptome

mkdir -p $ANALYSIS/transcriptome_bam
for fa in $ANALYSIS/fastq/*.fastq; do
	name=$(basename $fa .fastq)
	bowtie2 --very-sensitive -x $ANALYSIS/bowties_indexes/transcriptome -U $fa --threads 5 | samtools sort | samtools view -b -F 20  > $ANALYSIS/transcriptome_bam/${name}.bam
	samtools index $ANALYSIS/transcriptome_bam/${name}.bam $ANALYSIS/transcriptome_bam/${name}.bam.bai
done

mkdir -p $ANALYSIS/whole_transcriptome_bam
for fa in $ANALYSIS/fastq/*.fastq; do
	name=$(basename $fa .fastq)
	bowtie2 --very-sensitive -x $ANALYSIS/bowties_indexes_whole_transcriptome/transcriptome -U $fa --threads 5 | samtools sort | samtools view -b -F 20 > $ANALYSIS/whole_transcriptome_bam/${name}.bam
	samtools index $ANALYSIS/whole_transcriptome_bam/${name}.bam $ANALYSIS/whole_transcriptome_bam/${name}.bam.bai
done

(for bam in $ANALYSIS/transcriptome_bam/*.bam; do
	name=$(basename $bam .bam)
	echo -ne "$name\t"
	samtools view -F 276 $bam | wc -l 
done)>$ANALYSIS/transcriptome_depth.txt
