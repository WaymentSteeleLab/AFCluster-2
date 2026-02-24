#!/bin/bash
echo $HOSTNAME
export HOME=$PWD
mkdir -p ./cache/colabfold/params/
export XDG_CACHE_HOME=./cache
export MPLCONFIGDIR=./cache

# Try to copy weights from staging, fallback to downloading if unavailable
if [ -d "./af2_weights" ] && [ "$(ls -A ./af2_weights 2>/dev/null)" ]; then
    echo "Attempting to copy weights from staging..."
    if cp ./af2_weights/* ./cache/colabfold/params/ 2>/dev/null && [ "$(ls -A ./cache/colabfold/params/ 2>/dev/null)" ]; then
        echo "Successfully copied weights from staging"
    else
        echo "Failed to copy from staging, downloading weights..."
        python -m colabfold.download
    fi
else
    echo "Staging directory not available, downloading weights..."
    python -m colabfold.download
fi

# move into repo
cd AFCluster-2

NAME=$1
SEQ=$2

echo ">" $NAME > ${NAME}.fasta
echo ${SEQ} >> ${NAME}.fasta

python afcluster.py --input ${NAME}.fasta --config configs/afcluster.yml

mkdir ../results_a3m
cp -r output/${NAME}/clusters/* ../results_a3m/

mkdir ../results_pdb
cp -r output/${NAME}/preds/* ../results_pdb/
