import glob
import os
import shutil
import pandas as pd
import numpy as np
import sqlite3 as sl
import argparse

parser = argparse.ArgumentParser(description='writing snp metrics to database per-chromosome for reclustering')
parser.add_argument('--chrom', type=str, default='nope', help='chromosome string')

args = parser.parse_args()

chrom = args.chrom


idat_path = '/data/CARD/PD/GP2/raw_genotypes/GP2_recluster/GS_output'
swarm_scripts_dir = f'/data/CARD/PD/GP2/swarm_scripts'

samples_list = list(set([idat.split('/')[-1].replace('_Grn.idat','').replace('_Red.idat','') for idat in glob.glob(f'{idat_path}/*/*.idat')]))
barcodes_list = list(set([x.split('_')[0] for x in samples_list]))

con = sl.connect(f'{idat_path}/gp2_snp_metrics_chr{chrom}.db')
# cursor = con.cursor()

missing_metrics_list = []

for sample in samples_list:
    code = sample.split('_')[0]
    snp_metrics_path = f'{idat_path}/{code}/snp_metrics_{sample}_chr{chrom}.csv'
    if os.path.isfile(snp_metrics_path):
        snp_metrics = pd.read_csv(snp_metrics_path)
        snp_metrics.loc[:,['Sample_ID','snpID', 'BAlleleFreq', 'LogRRatio', 'R', 'Theta', 'GenTrain_Score', 'GT']].to_sql(name='METRICS', if_exists='append', index=False, con=con)
    else:
        missing_metrics_list.append(sample)

# now with last metrics file, create snp table (every metrics file has same snps)
snp_metrics.loc[:,['snpID', 'chromosome', 'position', 'a1', 'a2', 'maf']].to_sql(name='SNPS', if_exists='append', index=False, con=con)
con.commit()
missing_metrics_list_out = list(set(missing_metrics_list))
missing_metrics_df = pd.DataFrame({'missing': missing_metrics_list_out})
missing_metrics_df.to_csv(f'{idat_path}/gp2_snp_metrics_chr{chrom}.MISSING', header=False, index=False)