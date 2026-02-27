import os
import sys
import glob
import yaml
import argparse
import subprocess
import numpy as np
import pandas as pd

from src.cluster import *
from src.utils.msa import *
from src.utils.seqs import *
from src.utils.mmseqs import *

def get_labels(args, df):
    if args.cluster_method == "dbscan":
        return cluster_DBSCAN(args, df)
    raise ValueError(f"Unknown clustering method")

def run_cluster(args, subfolder, input):

    with open(f"{subfolder}/{args.keyword}.log", "w") as f:
        IDs, seqs_incl_del = load_fasta(input)
        seqs = clean_seqs(seqs_incl_del)

        df = pd.DataFrame({'SequenceName': IDs, 'sequence': seqs, 'seq_incl_del': seqs_incl_del})
        query_ = df.iloc[:1]
        df = df.iloc[1:]
        L = len(df.sequence.iloc[0])
        df['frac_gaps'] = [x.count('-') / L for x in df['sequence']]

        df = df.loc[df.frac_gaps < float(args.gap_cutoff)]
        f.write(f"Filtered sequences by gap_cutoff={args.gap_cutoff}\n")
        
        df, clusters = get_labels(args, df)

        f.write(f"Found {len(clusters)} clusters using {args.cluster_method}\n")

        cluster_dir = os.path.join(subfolder, "clusters",)
        os.makedirs(cluster_dir, exist_ok=True)
        for clust in clusters:
            tmp = df.loc[df.dbscan_label == clust]
            out = pd.concat([query_, tmp], axis=0)
            outpath = os.path.join(cluster_dir, f"{args.keyword}_{clust:03d}.a3m")
            write_fasta(out.SequenceName.tolist(), out.seq_incl_del.tolist(), outfile=outpath)
            f.write(f"Wrote {outpath} (n={len(out)})\n")
    os.remove(f"{subfolder}/{args.keyword}.log")

def run_neighborcluster(args, subfolder, input):
    '''Downsample original MSA and of those, take N closest sequences to be the clusters.
    Each prediction is for different variant'''

    IDs, seqs = load_fasta(input)
    if args.downsample:
        L = len(seqs[0])
        #filter by frac gaps
        gap_filt_inds = [i for i,x in enumerate(seqs) if x.count('-') / L < float(args.gap_cutoff)]

        filt_msa_file = f"{subfolder}/filt_gaps.a3m"
        write_fasta([IDs[i] for i in gap_filt_inds], [seqs[i] for i in gap_filt_inds], outfile=filt_msa_file)
        IDs, seqs = downsample_msa(filt_msa_file)
        print(f'downsampled to {len(filtered_IDs)} variants')

    counter=0
    cluster_dir = os.path.join(subfolder, "clusters",)
    os.makedirs(cluster_dir, exist_ok=True)

    for ind, seq in list(zip(IDs, seqs)):
        closest_inds = get_closest_n_seqs(seq, seqs, n=args.num_neighbors)
        outpath = os.path.join(cluster_dir, f"{args.keyword}_{counter:03d}.a3m")

        inds_for_cluster = [IDs[x] for x in closest_inds]
        seqs_for_cluster = [seqs[x] for x in closest_inds]
        seqs_for_cluster = fix_neighborcluster_msas(seqs_for_cluster)

        write_fasta(inds_for_cluster, seqs_for_cluster, outfile=outpath)
        counter+=1

def main(args):
    if args.input is not None:
        ids, seqs = load_fasta(args.input); print(ids)
    
    elif args.msa is not None:
        ids = [os.path.basename(args.msa).split('.')[0]]
        seqs = [''] #dummy seq

    for id_, seq_ in zip(ids, seqs): 
        args.keyword += f'_{id_}'
        subfolder = os.path.join(args.outdir, id_)
        os.makedirs(subfolder, exist_ok=True)

        if args.msa is not None:
            msa_file = args.msa
            print(f'Using provided MSA: {msa_file}')
        else:
            msa_file = os.path.join(subfolder, f'{id_}.a3m')

        if not os.path.exists(msa_file): 
            print(f'MSA file not found, generating...')
            msa_seqs = run_mmseqs(seq_, args.tmpdir) # I think this is msa lines?
            with open(msa_file, "w") as a3m:
                a3m.write(msa_seqs[0])
        
        if args.cluster_method=='dbscan':
            print(f'Running DBSCAN clustering...')
            run_cluster(args, subfolder, msa_file)
        elif args.cluster_method=='neighbor':
            print(f'Running neighborcluster')
            run_neighborcluster(args, subfolder, msa_file)

        print(f'Running structure prediction...')
        pred_dir = os.path.join(subfolder, 'preds')
        os.makedirs(pred_dir, exist_ok=True)
        for i in range(args.num_seeds):
            for fil in sorted(glob.glob(f"{subfolder}/clusters/*.a3m")):
                print(fil)
                fil_name = os.path.splitext(os.path.basename(fil))[0]
                os.makedirs(f'{pred_dir}/{fil_name}/s{i}', exist_ok=True)
                print(fil_name)
                if os.path.exists(f'{pred_dir}/{fil_name}/s{i}/{fil_name}_0.done.txt'):
                    continue

                sp_command = ['colabfold_batch', 
                                    '--use-dropout', 
                                    '--num-recycle', f'{args.n_recycles}', 
                                    '--random-seed', f'{i}', 
                                    '--jobname-prefix', f'{fil_name}',
                                    f'{fil}', f'{pred_dir}/{fil_name}/s{i}']
                if args.amber_relax:
                    sp_command.extend(['--amber', '--use-gpu-relax'])                   

                subprocess.run(sp_command)

        if args.zip_outputs:

            protein_ID = subfolder.split('/')[-1]
            shutil.copy(args.config, f'{subfolder}/config.yml')

            if os.path.exists('out2_rep_seq.fasta'): #from NeighborCluster
                shutil.copy('out2_rep_seq.fasta', f'{subfolder}/rep_seqs.fasta')
                os.remove('out2_rep_seq.fasta')
            shutil.make_archive(f'{args.keyword}', 'zip', args.outdir )

if __name__ == "__main__":
    p = argparse.ArgumentParser()

    p.add_argument("--input", type=str, help="Input fasta")
    p.add_argument('--msa', type=str, default=None, help="Input MSA")
    p.add_argument('--config', type=str, default='configs/afcluster.yml', help='config file')
   
    args = p.parse_args()

   #must provide either input or msa
    if args.input is None and args.msa is None:
        p.error("Either input or msa must be provided")

    if args.input is not None and args.msa is not None:
        p.error("Only one of input or msa can be provided")
    
    with open(args.config, "r") as f:
        cfg = yaml.safe_load(f)
    
    for k, v in cfg.items():
        setattr(args, k, v)

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.tmpdir, exist_ok=True)

    np.random.default_rng(seed=args.random_seed)
    
    os.environ['PATH'] = args.path_vars['PATH']
    os.environ["XDG_CACHE_HOME"] = args.path_vars['XDG_CACHE_HOME']
    os.environ["MPLCONFIGDIR"] = args.path_vars['MPLCONFIGDIR']
    
    main(args)
