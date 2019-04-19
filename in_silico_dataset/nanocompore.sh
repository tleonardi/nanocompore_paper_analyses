#!/bin/bash
export BASEDIR="/hps/nobackup/enright/tom/in_silico_dataset"
export DATA="/hps/nobackup/enright/nanopore/analyses/nanocompore_paper_analyses/in_silico_dataset/data/simulated_datasets2"
export REFERENCE="/hps/nobackup/enright/nanopore/analyses/nanocompore_paper_analyses/in_silico_dataset/data/references/random_guided_weight.fa"
export OUT="$BASEDIR/analysis"
mkdir -p $OUT

source activate $BASEDIR/conda_env
mkdir -p $OUT/nanocompore/
mkdir -p $OUT/logs/
bgadd -L 9 /insilico

for ds in $(seq -w 2 145); do
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
		--comparison_methods GMM,KS \
		--sequence_context 2 \
		--sequence_context_weights harmonic \
		--logit \
		--allow_warnings \
		--nthreads 20 \
		--overwrite"
	bsub -g /insilico -q standard -eo $OUT/logs/${samp_name}.err -oo $OUT/logs/${samp_name}.out -M 16000 -n 10 -R 'rusage[mem=16000]' $COMMAND
done


