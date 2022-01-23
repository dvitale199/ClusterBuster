# imports
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import streamlit as st

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








#### MAY WANT TO GRIDSEARCH OVER COVARIANCE TYPE!!!!
# def gtype_gmm(snp_theta_r_df, n_components, snp):
    
#     X = snp_theta_r_df[[theta_col, r_col]].copy()

#     gmm = GaussianMixture(
#         n_components=n_components,
#         covariance_type="diag",
#         random_state = 10).fit(X)

#     labels = gmm.predict(X)

#     out_dict = {
#         'gmm': gmm,
#         'X': X,
#         'y_pred': labels
#     }

#     # X.loc[:,'predicted_label'] = labels
#     return out_dict






def flagged_variants_dropdown(snp_list):
    pass



def search_variant():
    pass

def variant_report():
    pass

def recluster_variant():
    pass

def export_cluster():
    pass

def export_final_report():
    pass





################## TESTING: to be removed #################
    # report.drop(columns=['Index', 'Address', 'Chr', 'Position', 'GenTrain Score', 'Frac A', 'Frac C', 'Frac G', 'Frac T'], inplace=True)

    # Chop out to GType, Score and Theta dataframes.
    # Name_df = report[['Name']].copy()

    # GType_df = report.loc[:, report.columns.str.endswith(".GType")]
    # Theta_df = report.loc[:, report.columns.str.endswith(".Theta")]
    # R_df = report.loc[:, report.columns.str.endswith(".R")]

    # Name_df['Name.GType'] = Name_df.loc[:,'Name'] + ".GType"
    # Name_df['Name.Theta'] = Name_df.loc[:,'Name'] + ".Theta"
    # Name_df['Name.R'] = Name_df.loc[:,'Name'] + ".R"

    # # Merge names plus data types.

    # Name_GType_df = pd.concat([Name_df, GType_df], axis=1)
    # Name_GType_df.rename(columns={"Name.GType": "idx"}, inplace=True)
    # Name_GType_df.set_index('idx', inplace=True)
    # Name_GType_df.drop(columns=['Name','Name.Theta','Name.R'], inplace=True)
    # Name_GType_df.columns = Name_GType_df.columns.str.rstrip(".GType")

    # Name_Theta_df = pd.concat([Name_df, Theta_df], axis=1)
    # Name_Theta_df.rename(columns={"Name.Theta": "idx"}, inplace=True)
    # Name_Theta_df.set_index('idx', inplace=True)
    # Name_Theta_df.drop(columns=['Name','Name.GType','Name.R'], inplace=True)
    # Name_Theta_df.columns = Name_Theta_df.columns.str.rstrip(".Theta")

    # Name_R_df = pd.concat([Name_df, R_df], axis=1)
    # Name_R_df.rename(columns={"Name.R": "idx"}, inplace=True)
    # Name_R_df.set_index('idx', inplace=True)
    # Name_R_df.drop(columns=['Name','Name.Theta','Name.GType'], inplace=True)
    # Name_R_df.columns = Name_R_df.columns.str.rstrip(".R")

    # #  Transpose the data frames and make the names of the variants plus the suffixes the columns.
    # GType_transposed_df = Name_GType_df.transpose()
    # Theta_transposed_df = Name_Theta_df.transpose()
    # R_transposed_df = Name_R_df.transpose()

    # # Smash everything together and get ready for plotting.
    # temp_df = GType_transposed_df.merge(Theta_transposed_df, left_index=True, right_index=True)
    # clusterbuster_df = temp_df.merge(R_transposed_df, left_index=True, right_index=True)
    # snps_list = list(report['Name'].unique())
    # flagged_snps_list = list()
