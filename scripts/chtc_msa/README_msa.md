# NeighborCluster on CHTC, from a pre-built MSA

Folds NC clusters on CHTC starting from an `.a3m` you supply, so no MSA is built
at runtime.

## Files

| File | Purpose |
|---|---|
| `nc_msa.sub` | submit description; queues over MSA files |
| `run_nc_msa_chtc.sh` | executable; stages weights, calls the clusterer, folds one shard |
| `nc_cluster.py` | the NeighborCluster step — writes one cluster `.a3m` per center |
| `stage_af2_weights.sh/.sub` | one-time AF2 weights staging |

## 1. Stage the AF2 weights (once)

```bash
condor_submit stage_af2_weights.sub
cat stage_weights.out          # the verification block prints the osdf:/// path
```

Put that path in `nc_msa.sub`.

## 2. Set up

```bash
cd ~/AFCluster-2/scripts/chtc_msa
mkdir -p condor_log results
cp /path/to/chtc_msa/name.a3m .        # your MSA, alongside the .sub
```


```bash
condor_submit nc_msa.sub
```

## 4. Collect

`transfer_output_remaps` gives each shard its own directory:

```
results/name.a3m_shard0_pdb/
results/name.a3m_shard1_pdb/
```

Merge:

```bash
mkdir -p merged_pdb
cp -r results/*_pdb/* merged_pdb/
```

Each shard also returns the full cluster set in `*_a3m/`, identical across shards
since clustering is deterministic. Keep one, discard the rest.

---

## Sharding

NC builds one cluster per sequence.

`msa_jobs.csv` and `smoketest.csv` ship as templates naming `my_msa.a3m` —
edit them to your own MSA filename. The shipped 4 shards × 3 seeds is
≈ 6 h/job, inside the ~12 h `short` limit.

Shard *i* folds clusters `i, i+N, i+2N, …` — every cluster folded once, no
overlap. To change it:

```bash
N=4; SEEDS=3; MSA=name.a3m
for i in $(seq 0 $((N-1))); do echo "$MSA,$i,$N,$SEEDS"; done > msa_jobs.csv
```

Shards are independent and re-runnable; completed predictions are skipped on
retry.

## Runner arguments

```
run_nc_msa_chtc.sh <MSA> <SHARD> <NSHARDS> <SEEDS> [MAX_CLUSTERS] [NUM_NEIGHBORS] [GAP_CUTOFF]
```

The last three default to the values in `configs/afneighborcluster.yml`.
`MAX_CLUSTERS` caps the number of cluster centers.

## Clustering without a GPU

`nc_cluster.py` runs standalone, so cluster counts and depth can be checked
before submitting:

```bash
python scripts/chtc_msa/nc_cluster.py MSA.a3m NAME \
    --max-clusters 300 --num-neighbors 30 --gap-cutoff 0.7
```

Run it from the repo root, on the submit node or in the container — it needs
`biopython`, `numpy` and `pyyaml` via `src.utils.msa`.

## When a job fails

`condor_log/*.err`, and the `.out`, which echoes hostname, GPU, the resolved MSA
path and its sequence count, plus `n_recycles` and the cluster total.
