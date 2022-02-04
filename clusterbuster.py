# imports
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import streamlit as st

import plotly.express as px
import plotly.graph_objects as go


# from __future__ import absolute_import

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

import plotly.graph_objects as go


def calculate_maf(gtype_df):
    '''
    in:
        gtype_df (snpid    gtype1  gtype2 ... )
    out:
            'maf_df': minor allele frequency df (snpid  maf)
        '''
    
    gtypes_map = {
        'AA':0,
        'AB':1,
        'BA':1,
        'BB':2,
        'NC':np.nan
        }
    gtypes = gtype_df.set_index('snpid').replace(gtypes_map).copy()
    
    # count only called genotypes
    N = gtypes.shape[1]-gtypes.isna().sum(axis=1)
    freq = pd.DataFrame({'freq': gtypes.sum(axis=1)/(2*N)})
    freq.loc[:,'maf'] = np.where(freq < 0.5, freq, 1-freq)
    maf_out = freq.drop(columns=['freq']).reset_index()

    return maf_out

# @st.cache()
def parse_report(report_in, flag_maf, flag_gentrain):

    ''' 
    in: 
        reportfile path
    out: 
        dict {
            report: variant     theta     r 
            flagged_vars: variant   maf     gt_score    flag
            }
            
    '''
    
    report = report_in.drop(columns=['Index', 'Address', 'Chr', 'Position', 'GenTrain Score', 'Frac A', 'Frac C', 'Frac G', 'Frac T'])

    # Chop out to GType, Score and Theta dataframes.
    Name_df = report[['Name']].copy()
    GType_df = report.loc[:, report.columns.str.endswith(".GType")]
    Theta_df = report.loc[:, report.columns.str.endswith(".Theta")]
    R_df = report.loc[:, report.columns.str.endswith(".R")]

    Name_df['Name.GType'] = Name_df.loc[:,'Name'] + ".GType"
    Name_df['Name.Theta'] = Name_df.loc[:,'Name'] + ".Theta"
    Name_df['Name.R'] = Name_df.loc[:,'Name'] + ".R"

    # Merge names plus data types.
    Name_GType_df = pd.concat([Name_df, GType_df], axis=1).rename(columns={"Name.GType":"idx"}).drop(columns=['Name', 'Name.Theta', 'Name.R'])
    Name_GType_df.columns = Name_GType_df.columns.str.replace(".GType", "", regex=False)
    Name_GType_df_final = Name_GType_df.set_index('idx')

    Name_Theta_df = pd.concat([Name_df, Theta_df], axis=1).rename(columns={"Name.Theta": "idx"}).set_index('idx').drop(columns=['Name', 'Name.GType', 'Name.R'])
    Name_Theta_df.columns = Name_Theta_df.columns.str.replace(".Theta", "", regex=False)

    Name_R_df = pd.concat([Name_df, R_df], axis=1).rename(columns={"Name.R":"idx"}).set_index('idx').drop(columns=['Name', 'Name.Theta', 'Name.GType'])
    Name_R_df.columns = Name_R_df.columns.str.replace(".R", "", regex=False)

    #  Transpose the data frames and make the names of the variants plus the suffixes the columns.
    GType_transposed_df = Name_GType_df_final.transpose()
    Theta_transposed_df = Name_Theta_df.transpose()
    R_transposed_df = Name_R_df.transpose()

    # Smash everything together and get ready for plotting.
    temp_df = GType_transposed_df.merge(Theta_transposed_df, left_index=True, right_index=True)
    clusterbuster_df = temp_df.merge(R_transposed_df, left_index=True, right_index=True)
    clusterbuster_out = clusterbuster_df.reset_index().rename(columns={'index':'IID'})
    
    # get gentrain scores
    
    gtrain_scores_df = report_in.loc[:,['Name', 'Chr', 'Position','GenTrain Score','Frac A', 'Frac C', 'Frac G', 'Frac T']]
    gtrain_scores_df.columns = ['snpid','Chr','Position','gentrain_score','Frac_A','Frac_C','Frac_G','Frac_T']
    # calculate maf
    gtype_df = Name_GType_df.copy()
    gtype_df.loc[:,'snpid'] = gtype_df.loc[:,'idx'].str.replace(".GType", "", regex=False)
    gtype_to_maf = gtype_df.drop(columns=['idx'])
    maf_scores_df = calculate_maf(gtype_to_maf)
    flag_df = maf_scores_df.merge(gtrain_scores_df, how='inner', on='snpid')
    flag_df.loc[:,'maf_flag'] = np.where(flag_df.maf<flag_maf, True, False)
    flag_df.loc[:,'gentrain_flag'] = np.where(flag_df.gentrain_score<flag_gentrain, True, False)

    missing_df = pd.DataFrame()

    out_dict = {
        'clusterbuster_df': clusterbuster_out,
        'flagged_snps': flag_df, 
    }

    return out_dict


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


def plot_clusters(df, x_col, y_col, gtype_col, snpid):

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
        width=700, height=497)

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
        color="gtype_out", 
        color_discrete_map=cmap,
        range_x=xlim, range_y=ylim,
        width=938, height=630,
        labels={"gtype_out": ""}
        )

    # fig.data[0].showlegend = False
    # fig.data[3].showlegend = False
    # fig.data[6].showlegend = False

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











# def cluster_buster_plot(df_cluster, df_recluster, x_col, y_col, gtype_col, snpid):

#     d3 = px.colors.qualitative.D3

#     cmap = {
#         'AA': d3[0],
#         'AB': d3[1],
#         'BA': d3[1],
#         'BB': d3[2],
#         'NC': d3[3]
#     }

#     # gtypes_list = (df[gtype_col].unique())
#     xmin, xmax = df[x_col].min(), df[x_col].max()
#     ymin, ymax = df[y_col].min(), df[y_col].max()

#     xlim = [xmin-.1, xmax+.1]
#     ylim = [ymin-.1, ymax+.1]

#     fig1 = px.scatter(
#         df_cluster, 
#         x=x_col, y=y_col, 
#         color=gtype_col, 
#         color_discrete_map=cmap, 
#         width=600, height=400)

    # fig1.update_xaxes(range=xlim)
    # fig1.update_yaxes(range=ylim)


# def plot_gmm(df, x_col, y_col, gtype_col, gmm, snpid, n_std=3):

#     d3 = px.colors.qualitative.D3

#     cmap = {
#         'AA': d3[0],
#         'AB': d3[1],
#         'BA': d3[1],
#         'BB': d3[2],
#         'NC': d3[3]
#     }

#     gtypes_list = (df[gtype_col].unique())

#     fig = go.Figure()

#     for gtype in gtypes_list:
#         df_ = df.loc[df[gtype_col]==gtype]
#         x = df_.loc[:,x_col]
#         y = df_.loc[:,y_col]
#         color = cmap[gtype]
#         fig.add_trace(
#             go.Scatter(
#                 x=x, y=y, 
#                 mode="markers",
#                 name=gtype,
#                 marker = {'color':color},
#                 hovertemplate="Genotype=%s<br>Theta=%%{x}<br>R=%%{y}<extra></extra>"% gtype
#                 ))

#         fig.update_layout(title_text=f'{snpid} Cluster Plot')
#         fig.update_layout(legend_title_text = 'Genotype')
#         fig.update_xaxes(title_text=x_col)
#         fig.update_yaxes(title_text=y_col, tickangle=-90)

        
#         for n_std in range(1, n_std):
#             fig.add_shape(type='path',
#             path=confidence_ellipse(x, y, n_std=1.96*n_std, size=100), 
#             line_color=color,
#             fillcolor=color,
#             opacity=0.2)
    
#     return fig




