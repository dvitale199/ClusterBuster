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

metrics = pd.read_csv(metrics_path)

for sample in metrics.Sample_ID.unique():

    s_metrics = metrics.loc[metrics['Sample_ID']==sample]
    out_path_n = f'{out_path}_{sample}_CB.csv'

    if os.path.isfile(out_path_n):
        s_metrics.to_csv(out_path_n, mode='a', header=False, index=False)
    else:
        s_metrics.to_csv(out_path_n, header=True, index=False)
                    