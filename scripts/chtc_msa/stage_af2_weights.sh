#!/bin/bash
#
# Download AF2 params once and put them in /staging
#
# Run once:  condor_submit stage_af2_weights.sub
#
set -euo pipefail

export HOME=$PWD
mkdir -p cache output
export XDG_CACHE_HOME=$PWD/cache
export MPLCONFIGDIR=$PWD/cache

PARAMS=$PWD/cache/colabfold/params

echo "downloading AF2 params ..."
python -m colabfold.download

echo
echo "downloaded into $PARAMS:"
ls -lh "$PARAMS" || { echo "ERROR: params dir missing"; exit 1; }

# CHTC staging layout is /staging/<first letter of user>/<user>/
USER_NAME="${1:-}"
if [ -z "$USER_NAME" ]; then
    echo "ERROR: no username given."
    echo "  This script takes the username as its first argument; the submit"
    echo "  file supplies it via  Arguments = \$ENV(USER)."
    echo "  Run directly as:  ./stage_af2_weights.sh <username>"
    exit 1
fi
INITIAL="${USER_NAME:0:1}"
DEST="/staging/${INITIAL}/${USER_NAME}/af2_weights"
echo "user=$USER_NAME  ->  dest=$DEST"

echo
echo "staging to $DEST"
if ! mkdir -p "$DEST" 2>/dev/null; then
    echo "ERROR: cannot create $DEST"
    echo "  Does /staging/${INITIAL}/${USER_NAME} exist and is it writable?"
    ls -ld "/staging/${INITIAL}" "/staging/${INITIAL}/${USER_NAME}" 2>&1 || true
    exit 1
fi

cp -rL "$PARAMS"/. "$DEST"/

# ColabFold re-downloads unless this marker is present next to the params
if [ ! -f "$DEST/download_finished.txt" ]; then
    echo "marker absent after copy; creating it"
    touch "$DEST/download_finished.txt"
fi

echo
echo "verifying $DEST"
NPZ=$(ls "$DEST"/*.npz 2>/dev/null | wc -l)
echo "  .npz files      : $NPZ"
echo "  completion mark : $([ -f "$DEST/download_finished.txt" ] && echo present || echo MISSING)"
echo "  total size      : $(du -sh "$DEST" | cut -f1)"

if [ "$NPZ" -lt 5 ]; then
    echo "ERROR: expected at least 5 .npz param files, found $NPZ"
    exit 1
fi

echo
echo "DONE — point transfer_input_files at:"
echo "  osdf:///chtc/staging/${INITIAL}/${USER_NAME}/af2_weights"
