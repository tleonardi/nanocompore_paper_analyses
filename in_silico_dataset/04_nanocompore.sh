#!/bin/bash
# This script uses bsub to submit jobs to an LSF cluster

export BASEDIR="$(git rev-parse --show-toplevel)/in_silico_dataset/"
export DATA="$BASEDIR/data/simulated_datasets"
export REFERENCE="$BASEDIR/data/references/random_guided_weight.fa"
export OUT="$BASEDIR/analysis"
mkdir -p $OUT

mkdir -p $OUT/nanocompore/
mkdir -p $OUT/logs/
bgadd -L 9 /insilico
for ds in $(seq -w 2 595); do
	ref=$DATA/dataset_0001
	samp_name="dataset_0$ds"
	samp=$DATA/$samp_name

	COMMAND="nanocompore sampcomp  \
		--file_list1 $ref/reads_1.tsv,$ref/reads_2.tsv \
		--file_list2 $samp/reads_1.tsv,$samp/reads_2.tsv \
		--label1 Control \
		--label2 Dataset \
		--fasta $REFERENCE \
		--outpath $OUT/nanocompore/$samp_name \
		--comparison_methods GMM,KS,MW,TT \
		--sequence_context 2 \
		--sequence_context_weights harmonic \
		--logit \
		--allow_warnings \
		--nthreads 20 \
		--overwrite"
	bsub -g /insilico -q standard -eo $OUT/logs/${samp_name}.err -oo $OUT/logs/${samp_name}.out -M 16000 -n 10 -R 'rusage[mem=16000]' $COMMAND
done

cp $DATA/index.tsv $OUT

