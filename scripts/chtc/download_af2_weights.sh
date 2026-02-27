#!/bin/bash
export HOME=$PWD
mkdir output
mkdir cache
export XDG_CACHE_HOME=./cache
export MPLCONFIGDIR=./cache

python -m colabfold.download

ls -lhat ./cache/colabfold/params/*

mv ./cache/colabfold/params/*.txt /staging/$(USER)/af2_weights
mv ./cache/colabfold/params/*.npz /staging/$(USER)/af2_weights
mv ./cache/colabfold/params/LICENSE /staging/$(USER)/af2_weights

echo DONE
