#!/bin/bash
export BASEDIR="$(git rev-parse --show-toplevel)/in_silico_dataset/"
export OUT="$BASEDIR/analysis/roc_data/"
export CONDA_ENV="$BASEDIR/../conda_environments/in_silico.yml"
mkdir -p $OUT/logs
bgadd -L 3 /insilico

for ds in $(seq -w 2 595); do
	samp_name="dataset_0$ds"
	COMMAND="conda activate $CONDA_ENV &&
		 python $BASEDIR/calc_roc.py $samp_name > $OUT/$samp_name.txt"
	bsub -g /insilico -q standard -eo $OUT/logs/${samp_name}.err -oo $OUT/logs/${samp_name}.out -M 8000 -n 1 -R 'rusage[mem=8000]' "$COMMAND"
done

cat $OUT/dataset_*.txt > $OUT/all_datasets.txt
