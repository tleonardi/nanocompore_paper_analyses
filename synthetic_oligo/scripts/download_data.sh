#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
CONDA_ENVS="$BASEDIR/conda_environments/"
ANALYSIS="${BASEDIR}/synthetic_oligo_timp/analysis"

mkdir -p $ANALYSIS/data
wget -O $ANALYSIS/data/3primem6aoligo.tar.gz https://s3.amazonaws.com/timpnanopore/oxford/171208_m6A/3primem6aoligo.tar.gz
