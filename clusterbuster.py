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
from imblearn.under_sampling import RandomUnderSampler



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

  # drop snps where gentrain score
  snp_metrics = snp_metrics.loc[~snp_metrics['GenTrain_Score'].isna()]

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

  print(f'training model')
  cv_ = GridSearchCV(estimator=clf, param_grid=param_grid, cv=5)
  cv_.fit(X_train, y_train)

  return cv_


def test_gt_clf(X_test, y_test, clf):
  print(f'testing model')
  pred = clf.predict(X_test)
  acc = accuracy_score(y_test, pred)

  out_dict = {
      'pred': pred,
      'accuracy':acc
      }

  return out_dict


def recluster_gts(df, clf, label_encoder, min_prob=0.8):

  X_pred = df.loc[:,['Theta','R']]
  y_pred = clf.predict_proba(X_pred)

  pred_labels = [i.argmax() for i in y_pred]
  max_probs = [probs.max() for probs in y_pred]  
  out_gts = label_encoder.inverse_transform(pred_labels)

  df.loc[:,'new_gt'] = out_gts
  df.loc[:,'new_gt_label'] = pred_labels
  df.loc[:, 'max_prob'] = max_probs
  df.loc[:, 'reclustered'] = np.where(df.max_prob >= min_prob, True, False)

  return df


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

