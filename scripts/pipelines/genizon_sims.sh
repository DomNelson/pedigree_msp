#!/bin/bash

set -e

timestamp() {
    date +"%Y-%m-%d_%H-%M-%S"
}

if [ ! $1 ]
then
    echo Usage: $0 outdir
    exit 1
fi
# OUTDIR=$1/$( timestamp )
OUTDIR=$1
mkdir -p $OUTDIR

BASEDIR=$( realpath ../../ )
SCRIPTDIR=$BASEDIR/scripts
PEDFILE=$BASEDIR/data/Luke/Genizon_4149gen_nov2019.txt

N_SAMPLES=4134
POPSIZE=5000
LENGTH=3e8
MU=1e-9
MODEL=dtwf

OUTFILE=$OUTDIR/${N_SAMPLES}samples_${POPSIZE}popsize_${LENGTH}len_${MU}mu_${MODEL}model_$(timestamp).h5


cd $SCRIPTDIR

python pedigree_simulate.py  \
    --pedfile $PEDFILE       \
    --num_samples $N_SAMPLES \
    --outfile $OUTFILE       \
    --length $LENGTH         \
    --mu $MU                 \
    --model $MODEL           \
    --popsize $POPSIZE
