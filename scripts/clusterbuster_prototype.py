import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np
import math
import argparse

import plotly.express as px
import plotly.graph_objects as go

from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import balanced_accuracy_score, accuracy_score


parser = argparse.ArgumentParser(description='run clusterbuster on input metrics')
parser.add_argument('--metrics', type=str, default='nope', help='metrics file in')
parser.add_argument('--train_metrics', type=str, default='nope', help='metrics for training set in')
parser.add_argument('--chrom', type=str, default='nope', help='chromosome string. must run per-chromosome at this time to train proper baseline model on that particular chromosome')
parser.add_argument('--train_n', type=int, default=500, help='number of snps used for training baseline model')
parser.add_argument('--min_gentrain', type=float, default=0.5, help='minimum gentrain score for training baseline model (between 0 and 1)')
parser.add_argument('--min_prob', type=float, default=0.8, help='minimum probability for GT prediction to be accepted')
parser.add_argument('--out', type=str, default='nope', help='reclustered df out')

args = parser.parse_args()


metrics = args.metrics
train_metrics = args.train_metrics
chrom = args.chrom
train_n = args.train_n
min_gentrain = args.min_gentrain
min_prob = args.min_prob
out = args.out


def munge_train_test(df, test_size=0.2, random_state=123):

    X, y = df.loc[:, ['Theta','R']], df.loc[:, 'GT_label']

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=123)

    out_dict = {
      'X_train': X_train,
      'X_test': X_test,
      'y_train': y_train,
      'y_test': y_test,
      }

    return out_dict


def fit_gt_clf(X_train, y_train, clf=RandomForestClassifier(), param_grid={'criterion':['gini', 'entropy']}, cv=5):

    cv_ = GridSearchCV(estimator=clf, param_grid=param_grid, cv=cv)
    cv_.fit(X_train, y_train)

    return cv_


def test_gt_clf(X_test, y_test, clf):

    pred = clf.predict(X_test)
    acc = accuracy_score(y_test, pred)

    out_dict = {
      'pred': pred,
      'accuracy':acc
      }

    return out_dict


def predict_gts(df, clf, label_encoder, min_prob=0.8):
    df = df.copy()
    X_pred = df.loc[:,['Theta','R']]
    y_pred = clf.predict_proba(X_pred)
    y_pred_list = [l.tolist() for l in y_pred]
    pred_labels = [i.argmax() for i in y_pred]
    max_probs = [probs.max() for probs in y_pred]  
    out_gts = label_encoder.inverse_transform(pred_labels)

    df.loc[:,'new_gt'] = out_gts
    df.loc[:,'new_gt_label'] = pred_labels
    df.loc[:,'probs'] = y_pred_list
    df.loc[:, 'max_prob'] = max_probs
    df.loc[:, 'reclustered'] = np.where(((df['max_prob'] >= min_prob) & (df['R'] > 0.2)), True, False)
    df.loc[:,'gt_out'] = np.where(df['reclustered']==True, df['new_gt'], df['GT'])

    return df


def plot_clusters(df, x_col='Theta', y_col='R', gtype_col='GT', title='snp plot'):

    d3 = px.colors.qualitative.D3

    cmap = {
        'AA': d3[0],
        'AB': d3[1],
        'BA': d3[1],
        'BB': d3[2],
        'NC': d3[3]
    }

    # gtypes_list = (df[gtype_col].unique())
    xmin, xmax = df[x_col].min(), df[x_col].max()
    ymin, ymax = df[y_col].min(), df[y_col].max()

    xlim = [xmin-.1, xmax+.1]
    ylim = [ymin-.1, ymax+.1]

    fig = px.scatter(
        df, 
        x=x_col, y=y_col, 
        color=gtype_col, 
        color_discrete_map=cmap, 
        width=650, height=497)

    fig.update_xaxes(range=xlim, nticks=10)
    fig.update_yaxes(range=ylim, nticks=10)
    
    fig.update_layout(
        margin=dict(l=0, r=76, t=63, b=75),
    )
    fig.update_layout(legend=dict(
    orientation="h",
    yanchor="bottom",
    y=1,
    xanchor="right",
    x=1
))
    fig.update_layout(legend_title_text='')
    out_dict = {
        'fig': fig,
        'xlim': xlim,
        'ylim': ylim
    }
    
    fig.update_layout(title_text=title)
    

    return out_dict


def plot_hist_contour(df, x_col, y_col, gtype_col, xlim, ylim):

    d3 = px.colors.qualitative.D3

    cmap = {
        'AA': d3[0],
        'AB': d3[1],
        'BA': d3[1],
        'BB': d3[2],
        'NC': d3[3]
    }
    
    n_trace = len(df[gtype_col].unique())

    fig = px.density_contour(
        df,
        x=x_col, y=y_col, 
        marginal_x="histogram", 
        marginal_y="histogram",
        color=gtype_col, 
        color_discrete_map=cmap,
        range_x=xlim, range_y=ylim,
        width=900, height=630,
        labels={gtype_col: ""}
        )

    fig2 = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap, range_x=xlim, range_y=ylim, width=938, height=625)

    for i in range(n_trace):
        fig.add_trace(fig2.data[i])
    
    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=.75,
        xanchor="right",
        x=1
    ))

    return fig


snp_metrics_in = pd.read_csv(f'{metrics}')
snp_metrics = snp_metrics_in.loc[(~snp_metrics_in['GenTrain_Score'].isna()) & (~snp_metrics_in['Theta'].isna()) & (~snp_metrics_in['R'].isna())]
snp_metrics['maf'] = snp_metrics['maf'].fillna(0.0)
snp_metrics_train = pd.read_csv(f'{train_metrics}')

# now check if there are enough training snps, if not, use supplied metrics to grab training snps
if snp_metrics.loc[(snp_metrics['GenTrain_Score'] >= min_gentrain), 'GT'].value_counts()['AB'] < train_n:
    train_snps = snp_metrics_train.loc[(snp_metrics_train['GenTrain_Score'] >= min_gentrain) & (~snp_metrics_train['Theta'].isna()) & (~snp_metrics_train['R'].isna())]
else:
    train_snps = snp_metrics.loc[(snp_metrics['GenTrain_Score'] >= min_gentrain)]

# train baseline models: one for autosomes + X and one for Y and M
le = LabelEncoder()

if chrom in [str(i) for i in range(1,23)] + ['X']:
    aa_sample = train_snps.loc[train_snps.GT=='AA'].sample(n=train_n)
    ab_sample = train_snps.loc[train_snps.GT=='AB'].sample(n=train_n)
    bb_sample = train_snps.loc[train_snps.GT=='BB'].sample(n=train_n)

    baseline_train = pd.concat([aa_sample, ab_sample, bb_sample], ignore_index=True)

if chrom in ['Y','M']:
    aa_sample = train_snps.loc[train_snps.GT=='AA'].sample(n=train_n)
    bb_sample = train_snps.loc[train_snps.GT=='BB'].sample(n=train_n)
    
    baseline_train = pd.concat([aa_sample, bb_sample], ignore_index=True)

# transform labels for main training set (AA, AB, BB)
baseline_train.loc[:, 'GT_label'] = le.fit_transform(baseline_train.loc[:, 'GT'])

# train baseline model to use if all calls in NC snp are NC
baseline_train_test = munge_train_test(baseline_train, test_size=0.2, random_state=123)
baseline_clf = fit_gt_clf(X_train=baseline_train_test['X_train'], y_train=baseline_train_test['y_train'])
baseline_clf_test = test_gt_clf(X_test=baseline_train_test['X_test'], y_test=baseline_train_test['y_test'], clf=baseline_clf)

# get list of nc snps
nc_snps = list(snp_metrics.loc[snp_metrics['GT'] == 'NC', 'snpID'].unique())

total_recluster_df = pd.DataFrame()
count = 0
for selected_snp in nc_snps:
    count += 1
    cb_metrics = snp_metrics.loc[snp_metrics['snpID'] == selected_snp]
    
    if chrom in [str(i) for i in range(1,23)] + ['X']:
        gts = ['AA','AB','BB']
    if chrom in ['Y','M']:
        gts = ['AA','BB']

    # random oversample missing- CREATE FUNCTION FROM THIS
    cb_metrics_called = cb_metrics.loc[cb_metrics['GT'] != 'NC']
    gt_unique = list(cb_metrics_called['GT'].unique())

    gt_counts = cb_metrics_called['GT'].value_counts(sort=True)
    gt_count_max = gt_counts.max()
    gt_counts_df = pd.DataFrame({'GT': list(gt_counts.index),
                                'count': list(gt_counts)})


    missing_gts = [gt for gt in gts if gt not in gt_unique]
    missing_gt_counts_df = pd.DataFrame({'GT': missing_gts, 'count':[0]*len(missing_gts)})

    gt_counts_df = pd.concat([gt_counts_df, missing_gt_counts_df])
    gt_counts_df.loc[:, 'new_needed'] = gt_count_max - gt_counts_df.loc[:, 'count']

    sup_snps = pd.DataFrame()

    # if we have representation of all 3 GTs, just use those for training
    if len(gt_counts) == 3:
        cb_metrics_final = cb_metrics

    # supplement each class with some high-quality snps from training set
    # supplementary gts added based on chrom above
    if (len(gt_counts) > 0) & (len(gt_counts) <= 3):
        for gt in gts:

            sup_tmp = baseline_train.loc[baseline_train['GT']==gt]
            sup_snps = pd.concat([sup_snps, sup_tmp], ignore_index=True)

        # add these supplementary snps to cb_metrics for training snp-specific model
        cb_metrics_final = pd.concat([cb_metrics, sup_snps])


        training_snps_total = cb_metrics_final.loc[cb_metrics_final['GT']!='NC'].copy()

        # create labels
        le = LabelEncoder()
        training_snps_total.loc[:, 'GT_label'] = le.fit_transform(training_snps_total.loc[:, 'GT'])

        # train model for snp
        model_split = munge_train_test(training_snps_total, test_size=0.2, random_state=123)
        clf = fit_gt_clf(X_train=model_split['X_train'], y_train=model_split['y_train'])
        clf_test = test_gt_clf(X_test=model_split['X_test'], y_test=model_split['y_test'], clf=clf)

        # predict NCs
        pred_df = cb_metrics_final.loc[cb_metrics_final['GT']=='NC'].copy()
        recluster_snp_df = predict_gts(df=pred_df, clf=clf, label_encoder=le, min_prob=min_prob)

    else:
        pred_df = cb_metrics.loc[cb_metrics['GT']=='NC']
        # recluster no calls in autosomes and X chrom
        recluster_snp_df = predict_gts(df=pred_df, clf=baseline_clf, label_encoder=le, min_prob=min_prob)
    
    total_recluster_df = pd.concat([total_recluster_df, recluster_snp_df], ignore_index=True)

out_cols = ['chromosome','position','snpID','Sample_ID','R','Theta','GenTrain_Score','a1','a2','GT','maf','new_gt','probs','max_prob','reclustered','gt_out']
    
total_recluster_df[out_cols].to_csv(f'{out}', index=False)