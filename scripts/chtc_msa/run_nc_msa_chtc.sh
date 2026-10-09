#!/bin/bash
#
# NeighborCluster from a PRE-BUILT MSA, with option for sharding.
#
# Arguments: <MSA_FILE> <SHARD_INDEX> <N_SHARDS> <NUM_SEEDS> [MAX_CLUSTERS] [NUM_NEIGHBORS]
#   MSA_FILE      .a3m, transferred w/ job
#   SHARD_INDEX   0-based index of this shard
#   N_SHARDS      total number of shards (use 1 to fold everything in one job)
#   NUM_SEEDS     ColabFold seeds per cluster
#   MAX_CLUSTERS  how many cluster centers to use. Default is 1 per cluster.
#   NUM_NEIGHBORS cluster size. Generally 100-500.
#   GAP_CUTOFF    max gap fraction to keep a sequence (default: YAML value). Around 0.1.
#
set -euo pipefail

echo "host: $HOSTNAME"
# nvidia-smi is not on PATH in the ColabFold image; this is informational only.
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
else
    echo "nvidia-smi not on PATH (normal in this container) - GPU presence is"
    echo "reported by ColabFold itself further down as 'Running on GPU'"
fi

export HOME=$PWD
mkdir -p ./cache/colabfold/params/
export XDG_CACHE_HOME=$PWD/cache
export MPLCONFIGDIR=$PWD/cache

# Resolve this script's own directory BEFORE any cd, so the helper module can
# be found afterwards. On CHTC the files are transferred flat into the job
# scratch dir, so nc_cluster.py may sit next to this script or in $PWD.
# HTCondor copies the Executable to the scratch root but transfers the repo as
# a subdir, so the module usually arrives at AFCluster-2/scripts/chtc_msa/.
SCRIPT_DIR=""
for d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" \
         "$PWD" \
         "$PWD/AFCluster-2/scripts/chtc_msa"; do
    if [ -f "$d/nc_cluster.py" ]; then SCRIPT_DIR=$d; break; fi
done
if [ -z "$SCRIPT_DIR" ]; then
    echo "ERROR: nc_cluster.py not found. Looked beside the script, in $PWD," >&2
    echo "       and in $PWD/AFCluster-2/scripts/chtc_msa. Scratch contains:" >&2
    ls -la . >&2
    exit 1
fi
echo "nc_cluster.py: $SCRIPT_DIR"

MSA_FILE=$1
SHARD=${2:-0}
NSHARDS=${3:-1}
NUM_SEEDS=${4:-3}
MAX_CLUSTERS=${5:-0}
NUM_NEIGHBORS=${6:-0}
GAP_CUTOFF=${7:-0}

NAME=$(basename "$MSA_FILE" .a3m)
# Prevent name collisions from transfer_output_remaps.
if [ "$MAX_CLUSTERS" != "0" ] || [ "$NUM_NEIGHBORS" != "0" ] || [ "$GAP_CUTOFF" != "0" ]; then
    NAME="${NAME}_nb${NUM_NEIGHBORS}_gc${GAP_CUTOFF}_mc${MAX_CLUSTERS}"
fi
echo "msa=$MSA_FILE name=$NAME shard=$SHARD/$NSHARDS seeds=$NUM_SEEDS max_clusters=$MAX_CLUSTERS num_neighbors=$NUM_NEIGHBORS gap_cutoff=$GAP_CUTOFF"

# Create output dirs, supports better error surfacing
mkdir -p "$PWD/results_a3m" "$PWD/results_pdb"
RESULTS_A3M="$PWD/results_a3m"
RESULTS_PDB="$PWD/results_pdb"

# Echo the failing line number on any unexpected exit
trap 'echo "ERROR: exited $? at line $LINENO" >&2' ERR

# ---- AF2 weights ------------------------------------------------------------
# The previous job printed "Successfully copied weights from staging" and STILL
# downloaded ~4 GB: the old test only asked "is the target dir non-empty", which
# passes even when the params or ColabFold's completion marker are missing.
#
# ColabFold 1.5.5 re-downloads unless BOTH are present in params/:
# osdf:/// can't transfer directories, so you have to stage a tarball for this script to untar:
#   cd /staging/a/$USER && tar -cf af2_weights.tar -C af2_weights .
PARAMS=$PWD/cache/colabfold/params
mkdir -p "$PARAMS"

if [ -f ./af2_weights.tar ]; then
    echo "untarring af2_weights.tar -> $PARAMS"
    tar -xf ./af2_weights.tar -C "$PARAMS"
elif [ -d ./af2_weights ] && ls ./af2_weights/*.npz >/dev/null 2>&1; then
    # directory form (e.g. shared filesystem rather than osdf)
    echo "linking ./af2_weights -> $PARAMS"
    rm -rf "$PARAMS"; ln -s "$(cd ./af2_weights && pwd)" "$PARAMS"
elif [ -n "${AF2_STAGING:-}" ] && ls "$AF2_STAGING"/*.npz >/dev/null 2>&1; then
    echo "linking $AF2_STAGING -> $PARAMS"
    rm -rf "$PARAMS"; ln -s "$AF2_STAGING" "$PARAMS"
else
    echo "no staged params found. scratch contains:"; ls -la . | head -15
fi

NPZ=$(ls "$PARAMS"/*.npz 2>/dev/null | wc -l | tr -d " ")
if [ -f "$PARAMS/download_finished.txt" ]; then MARK=yes; else MARK=no; fi
echo "params check: npz=$NPZ marker=$MARK"

if [ "$MARK" = "yes" ] && [ "$NPZ" -ge 5 ]; then
    echo "weights OK"
else
    echo "weights incomplete; downloading"
    rm -rf "$PARAMS"; mkdir -p "$PARAMS"
    python -m colabfold.download
fi

# ---- locate the MSA ---------------------------------------------------------
# condor lands transferred files in the job scratch dir ($PWD), but the repo is
# a subdir, so resolve to an absolute path BEFORE cd'ing into it.
if [ -f "$PWD/$MSA_FILE" ]; then
    MSA_ABS="$PWD/$MSA_FILE"
elif [ -f "$MSA_FILE" ]; then
    MSA_ABS="$(cd "$(dirname "$MSA_FILE")" && pwd)/$(basename "$MSA_FILE")"
else
    echo "ERROR: MSA $MSA_FILE not found in $PWD"; ls -la; exit 1
fi
echo "resolved msa: $MSA_ABS  ($(grep -c '^>' "$MSA_ABS") sequences)"

cd AFCluster-2

# ---- 1. cluster only (no folding) -------------------------------------------
# num_seeds: 0 makes afcluster.py's prediction loop a no-op, so this call just
# writes the cluster .a3m files. Folding is driven explicitly below so we can
# shard it.
python "$SCRIPT_DIR/nc_cluster.py" "$MSA_ABS" "$NAME" \
    --max-clusters  "$MAX_CLUSTERS" \
    --num-neighbors "$NUM_NEIGHBORS" \
    --gap-cutoff    "$GAP_CUTOFF"

N_RECYCLES=$(python -c "import yaml; print(yaml.safe_load(open('configs/afneighborcluster.yml')).get('n_recycles', 3))")
echo "n_recycles=$N_RECYCLES"

CLUSTER_DIR=output/${NAME}/clusters
PRED_DIR=output/${NAME}/preds
mkdir -p "$PRED_DIR"

# ---- 2. fold this shard -----------------------------------------------------
# mapfile requires bash 4+; use a portable read loop instead. Also tolerate the
# no-match case rather than dying under `set -e`.
ALL=()
while IFS= read -r line; do
    ALL+=("$line")
done < <(find "$CLUSTER_DIR" -maxdepth 1 -name '*.a3m' | sort)
TOTAL=${#ALL[@]}
echo "total clusters: $TOTAL"
if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: no cluster .a3m files were written to $CLUSTER_DIR"
    ls -la "$CLUSTER_DIR" 2>&1 || true
    exit 1
fi

for ((i = SHARD; i < TOTAL; i += NSHARDS)); do
    fil="${ALL[$i]}"
    base=$(basename "$fil" .a3m)
    for ((s = 0; s < NUM_SEEDS; s++)); do
        out="$PRED_DIR/$base/s$s"
        if compgen -G "$out/*.pdb" > /dev/null; then
            echo "skip $base s$s (done)"; continue
        fi
        mkdir -p "$out"
        echo "folding $base seed $s"
        colabfold_batch --use-dropout \
            --num-recycle "$N_RECYCLES" \
            --random-seed "$s" \
            --jobname-prefix "$base" \
            "$fil" "$out"
    done
done

# ---- 3. collect -------------------------------------------------------------
# Only the clusters this shard folded are returned, so shards never collide.
cp -r "$CLUSTER_DIR"/* "$RESULTS_A3M"/ 2>/dev/null || true
if [ -d "$PRED_DIR" ]; then
    cp -r "$PRED_DIR"/* "$RESULTS_PDB"/ 2>/dev/null || true
fi
echo "done: $(ls "$RESULTS_PDB" 2>/dev/null | wc -l) prediction dirs"
