import pandas as pd
import argparse
import numpy as np
import shutil
import sys
import subprocess

parser = argparse.ArgumentParser(description='run clusterbuster on input metrics')
parser.add_argument('--metrics', type=str, default='nope', help='metrics file in')
parser.add_argument('--bfile', type=str, default='nope', help='prefix of plink bfile for individual sample (everything before .bed/.bim/.fam)')
parser.add_argument('--ped', type=str, default='nope', help='prefix of pre-cb ped file for individual sample (everything before .ped/.map)')
parser.add_argument('--out', type=str, default='nope', help='prefix of reclustered ped out (everything before .ped/.map)')

args = parser.parse_args()

# import this from genotools later on
def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))
    

metrics = args.metrics
bfile_in = args.bfile
ped_in = args.ped
out = args.out

fam_path = f'{bfile_in}.fam'
bim_path = f'{bfile_in}.bim'
map_path =  f'{ped_in}.map'
ped_path = f'{ped_in}.ped'

metrics_df = pd.read_csv(metrics)
sample_recluster = metrics_df.loc[metrics_df['reclustered']==True]

snp_map = pd.read_csv(
    map_path, 
    sep='\s+', 
    header=None, 
    names=['chrom','snpID','pos','bp'], 
    dtype={'chrom':str}
    )

fam = pd.read_csv(fam_path, sep='\s+', header=None, names=['fid','iid','pat','mat','sex','pheno'], dtype={'fid':str, 'iid':str, 'pat':str, 'mat':str, 'sex':str, 'pheno':str})

bim = pd.read_csv(bim_path, sep='\s+', header=None, names=['chrom','snpID','bp','pos','a1','a2'])

# recluster = sample_recluster.merge(bim, how='left', on='snpID')
# for now, just first row because we are using 1 sample per-fam. may need to be adjusted for broader use
fam_list = fam.iloc[0].to_list()

# # set alleles for ped output- this can probably be simplified to not require new columns. deal with this later
# recluster.loc[recluster['new_gt']=='AA', 'call1_new'] = recluster.loc[recluster['new_gt']=='AA', 'bim_a1']
# recluster.loc[recluster['new_gt']=='AA', 'call2_new'] = recluster.loc[recluster['new_gt']=='AA', 'bim_a1']
# recluster.loc[recluster['new_gt']=='AB', 'call1_new'] = recluster.loc[recluster['new_gt']=='AB', 'bim_a1']
# recluster.loc[recluster['new_gt']=='AB', 'call2_new'] = recluster.loc[recluster['new_gt']=='AB', 'bim_a2']
# recluster.loc[recluster['new_gt']=='BB', 'call1_new'] = recluster.loc[recluster['new_gt']=='BB', 'bim_a2']
# recluster.loc[recluster['new_gt']=='BB', 'call2_new'] = recluster.loc[recluster['new_gt']=='BB', 'bim_a2']


sample_recluster.loc[sample_recluster['new_gt']=='AA', 'call1_new'] = sample_recluster.loc[sample_recluster['new_gt']=='AA', 'a1']
sample_recluster.loc[sample_recluster['new_gt']=='AA', 'call2_new'] = sample_recluster.loc[sample_recluster['new_gt']=='AA', 'a1']
sample_recluster.loc[sample_recluster['new_gt']=='AB', 'call1_new'] = sample_recluster.loc[sample_recluster['new_gt']=='AB', 'a1']
sample_recluster.loc[sample_recluster['new_gt']=='AB', 'call2_new'] = sample_recluster.loc[sample_recluster['new_gt']=='AB', 'a2']
sample_recluster.loc[sample_recluster['new_gt']=='BB', 'call1_new'] = sample_recluster.loc[sample_recluster['new_gt']=='BB', 'a2']
sample_recluster.loc[sample_recluster['new_gt']=='BB', 'call2_new'] = sample_recluster.loc[sample_recluster['new_gt']=='BB', 'a2']

# read lines of ped file (one line per sample)
with open(ped_path) as ped_file:
    lines = ped_file.readlines()
    lines = [line.rstrip().split(' ') for line in lines]
ped_file.close()

snps = lines[0][6:]
gt_call1_list = snps[::2] 
gt_call2_list = snps[1::2]

snp_map.loc[:,'call1'] = gt_call1_list
snp_map.loc[:,'call2'] = gt_call2_list

# gt_snp_map = snp_map.merge(recluster, how='left', on='snpID')
gt_snp_map = snp_map.merge(sample_recluster, how='left', on='snpID')

gt_snp_map.loc[:,'call1_out'] = np.where(gt_snp_map['call1_new'].isna(), gt_snp_map['call1'], gt_snp_map['call1_new'])
gt_snp_map.loc[:,'call2_out'] = np.where(gt_snp_map['call2_new'].isna(), gt_snp_map['call2'], gt_snp_map['call2_new'])

# build ped with alternating call{1,2}_out values
ped_out_list = [None]*(len(list(gt_snp_map['call1_out']))+len(list(gt_snp_map['call2_out'])))
ped_out_list[::2] = gt_snp_map['call1_out']
ped_out_list[1::2] = gt_snp_map['call2_out']

ped_final_list = fam_list + ped_out_list

with open(f'{out}.ped', 'w') as ped_out_file:
    for p in ped_final_list:
        ped_out_file.write(f'{p} ')

ped_out_file.close()

shutil.copy(map_path, f'{out}.map')

metrics_df['snpID'].to_csv(f'{out}.snps', header=False, index=False)


# use plink to make bed
make_bed_cmd = f'plink --file {out} --make-bed --out {out}'
# !module load plink/1.9; {make_bed_cmd}
shell_do(make_bed_cmd)

# copy reclustered bim for reference later
shutil.copy(f'{out}.bim', f'{out}_re_alleles.bim')
re_bim = pd.read_csv(f'{out}.bim', sep='\s+', header=None, names=['chrom','snpID','bp','pos','a1','a2'])

# now overwrite reclustered bim with original alleles
bim_merge = re_bim.merge(bim, how='left', on=['chrom','snpID'])
bim_merge[['chrom','snpID','bp_x','pos_x','a1_y', 'a2_y']].to_csv(f'{out}.bim', sep='\t', header=False, index=False)