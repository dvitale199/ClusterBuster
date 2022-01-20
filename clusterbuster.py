# imports
import pandas as pd
import numpy as np

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


def read_report(reportfile, flag_maf, flag_gencall):

    ''' 
    in: 
        reportfile path
    out: 
        dict {
            report: variant     theta     r 
            flagged_vars: variant   maf     gt_score    flag
            }
            
    '''
    report_in = pd.read_csv(reportfile, engine='c', dtype={'Chr':str, 'position':int})
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

    # get gentrain scores
    gtrain_scores_df = report_in.loc[:,['Name','GenTrain Score']]
    gtrain_scores_df.columns = ['snpid','gentrain_score']
    # calculate maf
    gtype_df = Name_GType_df.copy()
    gtype_df.loc[:,'snpid'] = gtype_df.loc[:,'idx'].str.replace(".GType", "", regex=False)
    gtype_to_maf = gtype_df.drop(columns=['idx'])
    maf_scores_df = calculate_maf(gtype_to_maf)
    flag_df = maf_scores_df.merge(gtrain_scores_df, how='inner', on='snpid')
    flag_df.loc[:,'maf_flag'] = np.where(flag_df.maf<flag_maf, True, False)
    flag_df.loc[:,'gencall_flag'] = np.where(flag_df.gentrain_score<flag_gencall, True, False)


    out_dict = {
        'clusterbuster_df': clusterbuster_df,
        'flagged_snps': flag_df, 
    }

    return out_dict





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
