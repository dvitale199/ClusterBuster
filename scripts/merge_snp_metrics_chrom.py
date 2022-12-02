import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--metrics_path', type=str, default='nope', help='input path to snp metrics directory')
parser.add_argument('--samples', type=str, default='nope', help='text file with single column of sampleids, no header')
parser.add_argument('--chrom', type=str, default='nope', help='chromosome string')
parser.add_argument('--out', type=str, default='nope', help='output path')

args = parser.parse_args()

metrics_path = args.metrics_path
samples = args.samples
chrom = args.chrom
out_path = args.out

samples_df = pd.read_csv(samples, header=None, names=['sampleid'])
samples_list = list(samples_df.sampleid)

out_chr_df = pd.DataFrame()

for sample in samples_list:
    code = sample.split('_')[0]

    mfile = f'{metrics_path}/{code}/snp_metrics_{sample}_chr{chrom}.csv'

    mfile_df = pd.read_csv(mfile)              
    out_chr_df = pd.concat([out_chr_df, mfile_df], ignore_index=True)

out_chr_df.to_csv(out_path, index=False)