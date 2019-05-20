#!/bin/bash
export BASEDIR="/hps/nobackup/enright/tom/in_silico_dataset"
export OUT="$BASEDIR/analysis/roc_data/"
mkdir -p $OUT/logs
bgadd -L 3 /insilico

for ds in $(seq -w 2 145 ); do
	samp_name="dataset_0$ds"
	COMMAND="source activate $BASEDIR/conda_env &&
		 python $BASEDIR/calc_roc.py $samp_name > $OUT/$samp_name.txt"
	bsub -g /insilico -q standard -eo $OUT/logs/${samp_name}.err -oo $OUT/logs/${samp_name}.out -M 8000 -n 1 -R 'rusage[mem=8000]' "$COMMAND"
done

cat $OUT/dataset_*.txt > $OUT/all_datasets.txt
