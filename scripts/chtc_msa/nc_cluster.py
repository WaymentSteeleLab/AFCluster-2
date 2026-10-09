#!/usr/bin/env python3
"""NeighborCluster: write one cluster .a3m per center from a single MSA.

For each center sequence, take its `num_neighbors` nearest sequences and write that neighborhood as its own .a3m

Run from the AFCluster-2 repo root (it imports `src.utils` and reads
`configs/afneighborcluster.yml`).

    python scripts/chtc_msa/nc_cluster.py MSA.a3m NAME \\
        [--max-clusters N] [--num-neighbors N] [--gap-cutoff F] [--outdir DIR]

"""
import argparse
import os
import shutil
import sys

import numpy as np
import yaml

sys.path.insert(0, os.getcwd())
from src.utils.msa import load_fasta, write_fasta, fix_neighborcluster_msas
from src.utils.seqs import encode_seqs


def to_aligned(seqs):
    """a3m -> aligned columns: drop lowercase insertions, keep match states."""
    return ["".join(c for c in s if c.isupper() or c == "-") for s in seqs]


def filter_by_gaps(IDs, seqs, gap_cutoff):
    """Drop sequences that are mostly gaps. The query (index 0) always stays."""
    aligned = to_aligned(seqs)
    L = len(aligned[0])
    keep = [0] + [i for i in range(1, len(seqs))
                  if aligned[i].count("-") / L < gap_cutoff]
    return [IDs[i] for i in keep], [seqs[i] for i in keep]


def pick_centers(n_seqs, max_clusters):
    """Cluster centers, spread evenly over MSA order. 0/None keeps them all."""
    if max_clusters and max_clusters < n_seqs:
        centers = [0] + [int(round(x)) for x in
                         np.linspace(1, n_seqs - 1, max_clusters - 1)]
        return sorted(set(centers))
    return list(range(n_seqs))


def neighborhood(onehot_msa, k, n_nb):
    """Indices of center k plus its n_nb nearest sequences, center first."""
    d = np.linalg.norm(onehot_msa - onehot_msa[k], axis=1)
    idx = [int(x) for x in np.argsort(d, kind="stable")[: n_nb + 1]]
    if k in idx:
        idx.remove(k)
    return [k] + idx[:n_nb]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("msa", help="input .a3m (absolute path is safest)")
    ap.add_argument("name", help="run name; clusters land in <outdir>/<name>/clusters")
    ap.add_argument("--max-clusters", type=int, default=0,
                    help="cap the number of cluster centers (0 = one per sequence)")
    ap.add_argument("--num-neighbors", type=int, default=0,
                    help="neighbors per cluster (0 = take from config)")
    ap.add_argument("--gap-cutoff", type=float, default=0.0,
                    help="drop sequences with more than this gap fraction "
                         "(0 = take from config)")
    ap.add_argument("--config", default="configs/afneighborcluster.yml")
    ap.add_argument("--outdir", default=None,
                    help="override the config's outdir")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))

    IDs, seqs = load_fasta(args.msa)
    print(f"loaded {len(IDs)} sequences")
    if not seqs:
        sys.exit(f"ERROR: no sequences read from {args.msa}")

    gap_cutoff = (args.gap_cutoff if args.gap_cutoff > 0
                  else float(cfg.get("gap_cutoff", 0.25)))
    IDs, seqs = filter_by_gaps(IDs, seqs, gap_cutoff)
    print(f"after gap_cutoff={gap_cutoff}: {len(IDs)} sequences")
    if len(seqs) < 2:
        sys.exit(f"ERROR: gap_cutoff={gap_cutoff} left {len(seqs)} sequences; "
                 "raise it (ColabFold MSAs have a median column gap near 0.5)")

    n_nb = (args.num_neighbors if args.num_neighbors > 0
            else int(cfg.get("num_neighbors", 30)))
    if n_nb >= len(seqs):
        print(f"WARNING: num_neighbors={n_nb} >= sequences={len(seqs)}; "
              "every cluster will be the whole MSA")
    print(f"num_neighbors={n_nb}")

    centers = pick_centers(len(seqs), args.max_clusters)
    if args.max_clusters and args.max_clusters < len(seqs):
        print(f"subsampled to {len(centers)} cluster centers (of {len(seqs)})")
    else:
        print(f"using all {len(centers)} sequences as cluster centers")

    outdir = os.path.join(args.outdir or cfg.get("outdir", "output"),
                          args.name, "clusters")
    # Clear stale clusters: on a re-run with fewer centers the old files would
    # remain and the fold loop would silently fold a mix of both runs.
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir, exist_ok=True)

    aligned = to_aligned(seqs)
    ohe = encode_seqs(aligned, max_len=max(len(s) for s in aligned))

    for n, k in enumerate(centers):
        idx = neighborhood(ohe, k, n_nb)
        write_fasta([IDs[i] for i in idx],
                    fix_neighborcluster_msas([seqs[i] for i in idx]),
                    outfile=os.path.join(outdir, f"{args.name}_{n:04d}.a3m"))
    print(f"wrote {len(centers)} clusters to {outdir}")


if __name__ == "__main__":
    main()
