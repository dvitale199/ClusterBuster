# imports
import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np
import streamlit as st
import math

import plotly.express as px
import plotly.graph_objects as go

# from __future__ import absolute_import

# from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import LabelEncoder, StandardScaler
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.pipeline import Pipeline
# from sklearn.model_selection import KFold
# from sklearn.linear_model import LinearRegression
# from sklearn.model_selection import GridSearchCV
# from sklearn.metrics import balanced_accuracy_score

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold
# from sklearn.mixture import GaussianMixture
# from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import balanced_accuracy_score, accuracy_score



def calculate_maf(gtype_df):
    
    gtypes_map = {
        'AA': 0,
        'AB': 1,
        'BA': 1,
        'BB': 2,
        'NC': np.nan
        }

    gtypes = gtype_df.pivot(index='snpID', columns='Sample_ID', values='GT').replace(gtypes_map)
    
    # count only called genotypes
    N = gtypes.shape[1]-gtypes.isna().sum(axis=1)
    freq = pd.DataFrame({'freq': gtypes.sum(axis=1)/(2*N)})
    freq.loc[:,'maf'] = np.where(freq < 0.5, freq, 1-freq)
    maf_out = freq.drop(columns=['freq']).reset_index()

    return maf_out


def clusterbuster_munge(metrics_file_path, out_path):

  snp_metrics = pd.read_csv(metrics_file_path)

  alt_split = snp_metrics.loc[:,'Alt'].str.split(',', expand=True)
  snp_metrics.loc[:,'Alt1'], snp_metrics.loc[:,'Alt2'] = alt_split.loc[:,0], alt_split.loc[:,1]

  snp_metrics.loc[(snp_metrics['GType']=='AA') & (snp_metrics['ALLELE_A']==1), 'GT'] = 'BB'
  snp_metrics.loc[(snp_metrics['GType']=='AA') & (snp_metrics['ALLELE_A']==0), 'GT'] = 'AA'
  snp_metrics.loc[(snp_metrics['GType']=='BB') & (snp_metrics['ALLELE_B']==1), 'GT'] = 'BB'
  snp_metrics.loc[(snp_metrics['GType']=='BB') & (snp_metrics['ALLELE_B']==0), 'GT'] = 'AA'
  
  snp_metrics.loc[(snp_metrics['GType']=='AB'), 'GT'] = 'AB'
  snp_metrics.loc[(snp_metrics['GType']=='NC'), 'GT'] = 'NC'
  snp_metrics.loc[:,'GT'] = snp_metrics.loc[:,'GT'].fillna('NC')

  # drop snps where gentrain score, theta, and r isna
  snp_metrics = snp_metrics.loc[(~snp_metrics['GenTrain_Score'].isna()) & (~snp_metrics['Theta'].isna()) & (~snp_metrics['R'].isna())]

  snp_metrics.loc[snp_metrics['ALLELE_A']==0, 'a1'] = snp_metrics.loc[snp_metrics['ALLELE_A']==0,'Ref']
  snp_metrics.loc[snp_metrics['ALLELE_A']==1, 'a1'] = snp_metrics.loc[snp_metrics['ALLELE_A']==1,'Alt1']
  snp_metrics.loc[snp_metrics['ALLELE_B']==0, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==0,'Ref']
  snp_metrics.loc[snp_metrics['ALLELE_B']==1, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==1,'Alt1']
  snp_metrics.loc[snp_metrics['ALLELE_B']==2, 'a2'] = snp_metrics.loc[snp_metrics['ALLELE_B']==2,'Alt2']

  # calculate maf for full 
  maf_df = calculate_maf(snp_metrics)
  snp_metrics_full = snp_metrics.merge(maf_df, how='left', on='snpID')

  snp_metrics_full.drop(columns=['Unnamed: 0'], inplace=True)

  # write individual chr cleaned snp metrics to file
  snp_metrics_full[[
      'Sample_ID',
      'snpID',
      'chromosome', 
      'position',
      'maf',
      'GenTrain_Score',
      'BAlleleFreq',
      'LogRRatio',
      'R',
      'Theta',
      'GT',
      'a1',
      'a2']].to_csv(out_path, index=False)

  return snp_metrics_full


def training_random_sample(df, chrom, n=250):

  if chrom in ['Y','M']:
    aa_sample = df.loc[df.GT=='AA'].sample(n=n)
    bb_sample = df.loc[df.GT=='BB'].sample(n=n)
    cb_metrics = pd.concat([aa_sample, bb_sample], ignore_index=True)
  else:
    aa_sample = df.loc[df.GT=='AA'].sample(n=n)
    ab_sample = df.loc[df.GT=='AB'].sample(n=n)
    bb_sample = df.loc[df.GT=='BB'].sample(n=n)
    cb_metrics = pd.concat([aa_sample, ab_sample, bb_sample], ignore_index=True)

  return cb_metrics


def gt_conf_intervals(df, col, conf=0.95):

  conf_z_dict = {
    .80 :	1.282,
    .85 :	1.440,
    .90 :	1.645,
    .95	: 1.960,
    .99	: 2.576,
    .995 : 2.807,
    .999 : 3.291,
    }

  stats = df.groupby(['snpID','GT'])[col].agg(['mean', 'count', 'std'])

  z = conf_z_dict[conf]

  high = []
  low = []

  for i in stats.index:
    m, c, s = stats.loc[i]
    high.append(m + z*s/math.sqrt(c))
    low.append(m - z*s/math.sqrt(c))
  conf_str = str(conf).replace('0.','')
  stats[f'ci{conf_str}_high'] = high
  stats[f'ci{conf_str}_low'] = low

  return stats



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


def gt_classify(df, min_gentrain=0.8, min_prob=0.8, train_n=250):

  '''
  train_n: number of snps per GT (AA, AB, BB) for baseline models to use for no-call snps
  '''

  df = df.loc[(~df['GenTrain_Score'].isna()) & (~df['Theta'].isna()) & (~df['R'].isna())].copy()

  recluster_snps = list(df.loc[df['GT']=='NC','snpID'].unique())

  # eventually move this outside of this function... causing issue with random sampling
  train_snps = df.loc[(df['GenTrain_Score']>=min_gentrain) & (df['GT']!='NC')].copy()
  total_recluster_df = pd.DataFrame()

  # train baseline models: one for autosomes + X and one for Y and M
  baseline_le = LabelEncoder()
  baseline_le_alt = LabelEncoder()
  # baseline_train = training_random_sample(train_snps, chrom=chrom, n=250)
  
  aa_sample = df.loc[df.GT=='AA'].sample(n=train_n)
  ab_sample = df.loc[df.GT=='AB'].sample(n=train_n)
  bb_sample = df.loc[df.GT=='BB'].sample(n=train_n)

  baseline_train = pd.concat([aa_sample, ab_sample, bb_sample], ignore_index=True)

  # create alternate training set (AA and BB gts ONLY) and transform labels
  baseline_train_alt = baseline_train.loc[baseline_train['chromosome'].isin(['Y','M'])]
  baseline_train_alt.loc[:, 'GT_label'] = baseline_le_alt.fit_transform(baseline_train_alt.loc[:, 'GT'])

  # transform labels for main training set (AA, AB, BB)
  baseline_train.loc[:, 'GT_label'] = baseline_le.fit_transform(baseline_train.loc[:, 'GT'])
  
  # train one model across all gts on all chroms
  baseline_train_test = munge_train_test(baseline_train, test_size=0.2, random_state=123)
  baseline_clf = fit_gt_clf(X_train=baseline_train_test['X_train'], y_train=baseline_train_test['y_train'])
  baseline_clf_test = test_gt_clf(X_test=baseline_train_test['X_test'], y_test=baseline_train_test['y_test'], clf=baseline_clf)

  # train a second model for just chrY/M
  baseline_alt_train_test = munge_train_test(baseline_train_alt, test_size=0.2, random_state=123)
  baseline_alt_clf = fit_gt_clf(X_train=baseline_alt_train_test['X_train'], y_train=baseline_alt_train_test['y_train'])
  baseline_alt_clf_test = test_gt_clf(X_test=baseline_alt_train_test['X_test'], y_test=baseline_alt_train_test['y_test'], clf=baseline_alt_clf)

  for selected_snp in recluster_snps:

    cb_metrics = df.loc[df.snpID==selected_snp]
    chrom = str(cb_metrics.chromosome.unique()[0])
    
    # use different possible gts for autosomes + X vs M and Y
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
    # on hold for now--- going to balance classes first and see performance
    # if len(gt_counts) == 3:
    #   cb_metrics_final = cb_metrics

    if (len(gt_counts) > 0) & (len(gt_counts) <= 3):
      for gt in gts:

        # eventually move random sample outside of this to be able to run with smaller number of recluster snps
        n_sample = int(gt_counts_df.loc[gt_counts_df['GT']==gt,'new_needed'].iloc[0]) + 1
        sup_tmp = train_snps.loc[train_snps['GT']==gt].sample(n_sample)
        sup_snps = pd.concat([sup_snps, sup_tmp], ignore_index=True)

      cb_metrics_final = pd.concat([cb_metrics, sup_snps])


      training_snps_total = cb_metrics_final.loc[cb_metrics_final['GT']!='NC'].copy()

      # create labels
      le = LabelEncoder()
      training_snps_total.loc[:, 'GT_label'] = le.fit_transform(training_snps_total.loc[:, 'GT'])

      # train model for snp
      train_test_split = munge_train_test(training_snps_total, test_size=0.2, random_state=123)
      clf = fit_gt_clf(X_train=train_test_split['X_train'], y_train=train_test_split['y_train'])
      clf_test = test_gt_clf(X_test=train_test_split['X_test'], y_test=train_test_split['y_test'], clf=clf)

      # predict NCs
      pred_df = cb_metrics_final.loc[cb_metrics_final['GT']=='NC'].copy()
      recluster_snp_df = predict_gts(df=pred_df, clf=clf, label_encoder=le, min_prob=min_prob)

    else:
      pred_df = cb_metrics.loc[cb_metrics['GT']=='NC']

      if chrom in [str(i) for i in range(1,23)] + ['X']:
        # recluster no calls in autosomes and X chrom
        recluster_snp_df = predict_gts(df=pred_df, clf=baseline_clf, label_encoder=baseline_le, min_prob=min_prob)
      
      if chrom in ['Y', 'X']:
        # recluster no calls in Y and M chroms
        reclustered_snp_df = predict_gts(pred_df, clf=baseline_alt_clf, label_encoder=baseline_le_alt, min_prob=min_prob)

    # else:
    #   for gt in gts:
    #     n_sample = len(cb_metrics.Sample_ID.unique())
    #     sup_tmp = train_snps.loc[train_snps['GT']==gt].sample(n_sample)
    #     sup_snps = pd.concat([sup_snps, sup_tmp], ignore_index=True)

    #   cb_metrics_final = pd.concat([cb_metrics, sup_snps])

    total_recluster_df = pd.concat([total_recluster_df, recluster_snp_df])

  return total_recluster_df




def view_table_slice(df_in, max_rows=20, **st_dataframe_kwargs):
    """Display a subset of a DataFrame or Numpy Array to speed up app renders.
    
    Parameters
    ----------
    df : DataFrame | ndarray
        The DataFrame or NumpyArray to render.
    max_rows : int
        The number of rows to display.
    st_dataframe_kwargs : Dict[Any, Any]
        Keyword arguments to the st.dataframe method.
    """
    df = df_in.drop(columns='index').copy()

    n_rows = len(df)
    if n_rows <= max_rows:
        st.write(df)
    else:
        start = st.slider('Start row', 0, n_rows - max_rows)
        end = start + max_rows
        df = df[start:end]

        st.table(df, **st_dataframe_kwargs)
        st.text('Displaying rows %i to %i of %i.' % (start, end - 1, n_rows))


def plot_clusters(df, x_col='Theta', y_col='R', gtype_col='GT'):

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


@st.cache()
def process_cnv_reports(BAF, LRR, BIM, sample_id):

    # BAF reduce and transpose.
    BAF_temp_reduced = BAF.loc[BAF['V1'] == sample_id]
    BAF_temp_reduced.drop(columns=['V1','V2','V3','V4','V5','V6'], inplace=True)
    BAF_transposed = BAF_temp_reduced.transpose()
    BAF_transposed.columns = ["BAF"]

    # LRR reduce and transpose.
    LRR_temp_reduced = LRR.loc[LRR['V1'] == sample_id]
    LRR_temp_reduced.drop(columns=['V1','V2','V3','V4','V5','V6'], inplace=True)
    LRR_transposed = LRR_temp_reduced.transpose()
    LRR_transposed.columns = ["LRR"]

    BAF_transposed.reset_index(drop=True, inplace=True)
    LRR_transposed.reset_index(drop=True, inplace=True)
    BIM.reset_index(drop=True, inplace=True)

    out_df = pd.concat([BAF_transposed, LRR_transposed, BIM], axis=1)

    return out_df


@st.cache
def csv_convert_df(df):
    
    return df.to_csv().encode('utf-8')



def gtype_gmm(snp_theta_r_df, n_components):
    IIDs = snp_theta_r_df.loc[:,'IID']
    X = snp_theta_r_df[['Theta','R']].copy()

    gmm = GaussianMixture(
        n_components=n_components,
        covariance_type="full",
        random_state = 10).fit(X)
    
    labels = gmm.predict(X)

    out_dict = {
        'gmm': gmm,
        'X': X,
        'y_pred': labels,
        'IID':IIDs
    }

    return out_dict


def parse_report(report_in, flag_maf, flag_gentrain):

    gtype_df = report_in
    maf_df = report_in.loc[:,['Sample_ID','snpID','GType']] 
    maf_scores_df = calculate_maf(maf_df)
    out_df = maf_scores_df.merge(gtype_df, how='inner', on='snpID')
    out_df.loc[:,'maf_flag'] = np.where(out_df.maf<flag_maf, True, False)
    out_df.loc[:,'gentrain_flag'] = np.where(out_df.GenTrain_Score<flag_gentrain, True, False)
    
    flagged_snps_list = list(out_df.loc[(out_df.gentrain_flag | out_df.maf_flag), 'snpID'].unique())

    # prune flagged_snps_list for those eligible for reclustering
    # snps that do not having missing theta or R 
    # snps with NC
    # snps where called genotypes have > 2 instances

    out_dict = {
      'cb_df': out_df,
      'flagged_snps': flagged_snps_list 
    }
    return out_df

def prep_theta_r_df(df, snpid):
  # Let's treat this as a simple regression problem
  df_ = df.loc[df.GType!='NC']
  pred_df = df.loc[df.GType=='NC']
  X_pred = pred_df.loc[:,['Theta','R']]

  X, y_ = df_.loc[:, ['Theta','R']], df_.loc[:,['GType']]

  le = LabelEncoder()
  y = le.fit_transform(y_.GType)

  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=123)
  
  out_dict = {
      'snpid': snpid,
      'X_train': X_train,
      'X_test': X_test,
      'y_train': y_train,
      'y_test': y_test,
      'X_pred': X_pred,
      'pred_df': pred_df,
      'labelencoder':le
  }
  
  return out_dict

    
def fit_snp_clf(X_train, X_test, y_train, y_test, snpid):

  clf = RandomForestClassifier()

  regressor = clf.fit(X_train, y_train)
  pred = regressor.predict(X_test)
  acc = balanced_accuracy_score(y_test,pred)
  
  out_dict = {
      'snpid':snpid,
      'regressor': regressor,
      'accuracy': acc
      }

  return out_dict


def predict_genotype_prob(X_pred, regressor, snpid, prob_cutoff=0.75):
    
  y_pred = regressor.predict_proba(X_pred) 
  
  out_dict = {
      'snp':snpid,
      'y_pred': y_pred
      }
  
  return out_dict


def cb_recluster_variant(df, snpid, min_prob=0.75):

  cb_data = prep_theta_r_df(df, snpid=snpid)
  X_train, X_test, y_train, y_test, X_pred = cb_data['X_train'], cb_data['X_test'], cb_data['y_train'], cb_data['y_test'], cb_data['X_pred']

  X_pred = X_pred.loc[~(X_pred['Theta'].isna() | X_pred['R'].isna())]
  
  gtype_clf = fit_snp_clf(X_train, X_test, y_train, y_test, snpid=snpid)
  regressor = gtype_clf['regressor']
  clf_acc = gtype_clf['accuracy']

  pred_gtype = predict_genotype_prob(X_pred=X_pred, regressor=regressor, snpid=snpid, prob_cutoff=0.75)
 
  pred_labels = [i.argmax() for i in pred_gtype['y_pred']]
  max_probs = [probs.max() for probs in pred_gtype['y_pred']]  
  out_gts = cb_data['labelencoder'].inverse_transform(pred_labels)

  pred_df = cb_data['pred_df'].copy()
  pred_df.loc[:,'new_gt'] = out_gts
  pred_df.loc[:,'new_gt_label'] = pred_labels
  pred_df.loc[:, 'max_prob'] = max_probs
  pred_df.loc[:, 'reclustered'] = np.where(pred_df.max_prob >= min_prob, True, False)
  pred_df.loc[:, 'clf_accuracy'] = clf_acc
  
  return pred_df

        
    #     for n_std in range(1, n_std):
    #         fig.add_shape(type='path',
    #         path=confidence_ellipse(x, y, n_std=1.96*n_std, size=100), 
    #         line_color=color,
    #         fillcolor=color,
    #         opacity=0.2)
    
    # return fig







######### old method to be removed 
# def parse_report(report_in, flag_maf, flag_gentrain):

#     ''' 
#     in: 
#         reportfile path
#     out: 
#         dict {
#             report: variant     theta     r 
#             flagged_vars: variant   maf     gt_score    flag
#             }
            
#     '''
    
#     report = report_in.drop(columns=['Index', 'Address', 'Chr', 'Position', 'GenTrain Score', 'Frac A', 'Frac C', 'Frac G', 'Frac T'])

#     # Chop out to GType, Score and Theta dataframes.
#     Name_df = report[['Name']].copy()
#     GType_df = report.loc[:, report.columns.str.endswith(".GType")]
#     Theta_df = report.loc[:, report.columns.str.endswith(".Theta")]
#     R_df = report.loc[:, report.columns.str.endswith(".R")]

#     Name_df['Name.GType'] = Name_df.loc[:,'Name'] + ".GType"
#     Name_df['Name.Theta'] = Name_df.loc[:,'Name'] + ".Theta"
#     Name_df['Name.R'] = Name_df.loc[:,'Name'] + ".R"

#     # Merge names plus data types.
#     Name_GType_df = pd.concat([Name_df, GType_df], axis=1).rename(columns={"Name.GType":"idx"}).drop(columns=['Name', 'Name.Theta', 'Name.R'])
#     Name_GType_df.columns = Name_GType_df.columns.str.replace(".GType", "", regex=False)
#     Name_GType_df_final = Name_GType_df.set_index('idx')

#     Name_Theta_df = pd.concat([Name_df, Theta_df], axis=1).rename(columns={"Name.Theta": "idx"}).set_index('idx').drop(columns=['Name', 'Name.GType', 'Name.R'])
#     Name_Theta_df.columns = Name_Theta_df.columns.str.replace(".Theta", "", regex=False)

#     Name_R_df = pd.concat([Name_df, R_df], axis=1).rename(columns={"Name.R":"idx"}).set_index('idx').drop(columns=['Name', 'Name.Theta', 'Name.GType'])
#     Name_R_df.columns = Name_R_df.columns.str.replace(".R", "", regex=False)

#     #  Transpose the data frames and make the names of the variants plus the suffixes the columns.
#     GType_transposed_df = Name_GType_df_final.transpose()
#     Theta_transposed_df = Name_Theta_df.transpose()
#     R_transposed_df = Name_R_df.transpose()

#     # Smash everything together and get ready for plotting.
#     temp_df = GType_transposed_df.merge(Theta_transposed_df, left_index=True, right_index=True)
#     clusterbuster_df = temp_df.merge(R_transposed_df, left_index=True, right_index=True)
#     clusterbuster_out = clusterbuster_df.reset_index().rename(columns={'index':'IID'})
    
#     # get gentrain scores
    
#     gtrain_scores_df = report_in.loc[:,['Name', 'Chr', 'Position','GenTrain Score','Frac A', 'Frac C', 'Frac G', 'Frac T']]
#     gtrain_scores_df.columns = ['snpid','Chr','Position','gentrain_score','Frac_A','Frac_C','Frac_G','Frac_T']
#     # calculate maf
#     gtype_df = Name_GType_df.copy()
#     gtype_df.loc[:,'snpid'] = gtype_df.loc[:,'idx'].str.replace(".GType", "", regex=False)
#     gtype_to_maf = gtype_df.drop(columns=['idx'])
#     maf_scores_df = calculate_maf(gtype_to_maf)
#     flag_df = maf_scores_df.merge(gtrain_scores_df, how='inner', on='snpid')
#     flag_df.loc[:,'maf_flag'] = np.where(flag_df.maf<flag_maf, True, False)
#     flag_df.loc[:,'gentrain_flag'] = np.where(flag_df.gentrain_score<flag_gentrain, True, False)

#     missing_df = pd.DataFrame()

#     out_dict = {
#         'clusterbuster_df': clusterbuster_out,
#         'flagged_snps': flag_df, 
#     }

#     return out_dict

