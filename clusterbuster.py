# imports
import pandas as pd
import numpy as np

def calculate_maf(gtype_df):
    '''
    in:
        gtype_df (id    gtype1  gtype2 ... )
    out:
            'maf_df': minor allele frequency df (id  maf)
        '''
    
    gtypes_map = {
        'AA':0,
        'AB':1,
        'BA':1,
        'BB':2,
        'NC':np.nan
        }
    gtypes = gtype_df.set_index('id').replace(gtypes_map).copy()
    
    # count only called genotypes
    N = gtypes.shape[1]-gtypes.isna().sum(axis=1)
    freq = pd.DataFrame({'freq': gtypes.sum(axis=1)/(2*N)})
    freq.loc[:,'maf'] = np.where(freq < 0.5, freq, 1-freq)
    maf_out = freq.drop(columns=['freq']).reset_index()

    return maf_out


def read_report(reportfile):

    ''' 
    in: 
        reportfile path
    out: 
        dict {
            report: variant     theta     r 
            flagged_vars: variant   maf     gt_score    flag
            }
            
    '''

    report = pd.read_csv(reportfile, engine='c', sep='\t')
    # report.drop(columns=['Index', 'Address', 'Chr', 'Position', 'GenTrain Score', 'Frac A', 'Frac C', 'Frac G', 'Frac T'], inplace=True)

    # Chop out to GType, Score and Theta dataframes.
    Name_df = report[['Name']].copy()

    GType_df = report.loc[:, report.columns.str.endswith(".GType")]
    Theta_df = report.loc[:, report.columns.str.endswith(".Theta")]
    R_df = report.loc[:, report.columns.str.endswith(".R")]

    Name_df['Name.GType'] = Name_df.loc[:,'Name'] + ".GType"
    Name_df['Name.Theta'] = Name_df.loc[:,'Name'] + ".Theta"
    Name_df['Name.R'] = Name_df.loc[:,'Name'] + ".R"

    # Merge names plus data types.

    Name_GType_df = pd.concat([Name_df, GType_df], axis=1)
    Name_GType_df.rename(columns={"Name.GType": "idx"}, inplace=True)
    Name_GType_df.set_index('idx', inplace=True)
    Name_GType_df.drop(columns=['Name','Name.Theta','Name.R'], inplace=True)
    Name_GType_df.columns = Name_GType_df.columns.str.rstrip(".GType")

    Name_Theta_df = pd.concat([Name_df, Theta_df], axis=1)
    Name_Theta_df.rename(columns={"Name.Theta": "idx"}, inplace=True)
    Name_Theta_df.set_index('idx', inplace=True)
    Name_Theta_df.drop(columns=['Name','Name.GType','Name.R'], inplace=True)
    Name_Theta_df.columns = Name_Theta_df.columns.str.rstrip(".Theta")

    Name_R_df = pd.concat([Name_df, R_df], axis=1)
    Name_R_df.rename(columns={"Name.R": "idx"}, inplace=True)
    Name_R_df.set_index('idx', inplace=True)
    Name_R_df.drop(columns=['Name','Name.Theta','Name.GType'], inplace=True)
    Name_R_df.columns = Name_R_df.columns.str.rstrip(".R")

    #  Transpose the data frames and make the names of the variants plus the suffixes the columns.
    GType_transposed_df = Name_GType_df.transpose()
    Theta_transposed_df = Name_Theta_df.transpose()
    R_transposed_df = Name_R_df.transpose()

    # Smash everything together and get ready for plotting.
    temp_df = GType_transposed_df.merge(Theta_transposed_df, left_index=True, right_index=True)
    clusterbuster_df = temp_df.merge(R_transposed_df, left_index=True, right_index=True)
    snps_list = list(report['Name'].unique())
    flagged_snps_list = list()

    out_dict = {

        'clusterbuster_df': clusterbuster_df,
        'flagged_snps': flagged_snps_list, 
        'all_snps': snps_list
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






