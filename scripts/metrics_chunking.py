import pandas as pd
import os
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='writing snp metrics to database per-chromosome for reclustering')
parser.add_argument('--metrics', type=str, default='nope', help='metrics file in')
parser.add_argument('--out', type=str, default='nope', help='metrics file out')

args = parser.parse_args()

metrics_path = args.metrics
out_path = args.out

out_cols = ['chromosome','position','snpID','Sample_ID','R','Theta','GenTrain_Score','a1','a2','GT','maf']

metrics = pd.read_csv(metrics_path)
n_snps = metrics.shape[0]
n_splits = n_snps/200
chrom_split = np.array_split(metrics, n_splits)

for n, split in enumerate(chrom_split):

    split_df = split.loc[:,out_cols]
    out_path_n = f'{out_path}_{n}.csv'

    if os.path.isfile(out_path_n):
        split_df.to_csv(out_path_n, mode='a', header=False, index=False)
    else:
        split_df.to_csv(out_path_n, header=True, index=False)
                    