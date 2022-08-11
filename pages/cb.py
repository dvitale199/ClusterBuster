from multiprocessing.connection import wait
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from sklearn.mixture import GaussianMixture
import time
from PIL import Image
from plotly.subplots import make_subplots
# from ClusterBuster.pages import gt
# from ClusterBuster.pages impo, gt
# from ClusterBuster import cnv
from pages import cb, cnv, gt

from clusterbuster import parse_report, view_table_slice, plot_clusters, gtype_gmm, plot_hist_contour, csv_convert_df

def run():
    st.sidebar.title('DTi Genotype Analysis')

    st.sidebar.subheader('Upload GenomeStudio Variant Report')

    upload_expander = st.sidebar.expander("Upload Report", expanded=False)

    with upload_expander:

        snp_metrics = st.file_uploader('Report')
        
        input_snp_file = st.file_uploader('Variant List (Optional)')

        if input_snp_file:
            input_snp_df = pd.read_csv(input_snp_file, header=None, names=['snpID'])
            input_snp_list = list(input_snp_df.loc[:,'snpID'])
        ##### Add try except to check file format
        #### add check to see if anything in the input list is not in snp data
        else:
            input_snp_list = None
        
    if snp_metrics:
        threshold_expander = st.sidebar.expander("Genotype Quality Control Thresholds (Advanced)", expanded=False)
        
        with threshold_expander:
            maf_threshold = st.slider('MAF Threshold', value=0.010, step=0.005, min_value=0.000, max_value=1.0)
            gentrain_threshold = st.slider('GenTrain Threshold', value=0.500, step=0.005, min_value=0.000, max_value=1.0)
        
        process_report_button = st.sidebar.button("Analyze Clusters")

        if process_report_button:
            report_in = pd.read_csv(snp_metrics, engine='c', dtype={'chromosome':str, 'position':int})
            report = parse_report(report_in, flag_maf=maf_threshold, flag_gentrain=gentrain_threshold)
            

            if report not in st.session_state:
                st.session_state['report'] = report
                st.session_state['update'] = False

            if report in st.session_state:
                st.session_state['report'] = report
                st.session_state['update'] = True

            with st.spinner('Calculating MAFs...'):
                time.sleep(2)
            with st.spinner('Flagging Variants for low MAF/GenTrain...'):
                time.sleep(2)

            st.sidebar.success('SNPs of interest have been flagged for low maf and gentrain')

    if 'report' in st.session_state.keys():
        
        report = st.session_state['report']
        cb_df = report
        snps_df = report['flagged_snps']
        snps_list = list(report.snpID.unique())




# keeping old method here for now
# def run():
    
#     st.sidebar.title('DTi Genotype Analysis')

#     st.sidebar.subheader('Upload GenomeStudio Variant Report')

#     upload_expander = st.sidebar.expander("Upload Report", expanded=False)

#     with upload_expander:

#         reportfile = st.file_uploader('Report')

#         input_snp_file = st.file_uploader('Variant List (Optional)')

#         if input_snp_file:
#             input_snp_df = pd.read_csv(input_snp_file, header=None, names=['snp'])
#             input_snp_list = list(input_snp_df.loc[:,'snp'])
#         ##### Add try except to check file format
#         #### add check to see if anything in the input list is not in snp data
#         else:
#             input_snp_list = None

#     if reportfile:
#         threshold_expander = st.sidebar.expander("Genotype Quality Control Thresholds (Advanced)", expanded=False)
        
#         with threshold_expander:
#             maf_threshold = st.slider('MAF Threshold', value=0.010, step=0.005, min_value=0.000, max_value=1.0)
#             gentrain_threshold = st.slider('GenTrain Threshold', value=0.500, step=0.005, min_value=0.000, max_value=1.0)
        
#         process_report_button = st.sidebar.button("Analyze Clusters")

#         if process_report_button:
#             report_in = pd.read_csv(reportfile, engine='c', dtype={'Chr':str, 'position':int})
#             report = parse_report(report_in, flag_maf=maf_threshold, flag_gentrain=gentrain_threshold)

#             if report not in st.session_state:
#                 st.session_state['report'] = report
#                 st.session_state['update'] = False

#             if report in st.session_state:
#                 st.session_state['report'] = report
#                 st.session_state['update'] = True

#             with st.spinner('Calculating MAFs...'):
#                 time.sleep(2)
#             with st.spinner('Flagging Variants for low MAF/GenTrain...'):
#                 time.sleep(2)

#             st.sidebar.success('SNPs of interest have been flagged for low maf and gentrain')






#     if 'report' in st.session_state.keys():
        
#         report = st.session_state['report']
#         cb_df = report['clusterbuster_df']
#         snps_df = report['flagged_snps']
#         snps_list = list(snps_df.snpid.unique())

#         # flagged dfs
#         total_flagged = snps_df[snps_df.maf_flag | snps_df.gentrain_flag].reset_index()
#         both_flagged = snps_df[snps_df.maf_flag & snps_df.gentrain_flag].reset_index()
#         maf_flagged = snps_df.loc[snps_df.maf_flag].reset_index()
#         gentrain_flagged = snps_df.loc[snps_df.gentrain_flag].reset_index()


#         snps_list = list(snps_df.snpid.unique())
#         flagged_snps_list = list(total_flagged.snpid.unique())
        
#         if total_flagged.shape[0]>0:
#             st.session_state['flagged'] = total_flagged
#         else:
#             total_flagged = None
#         # snp_list_upload_expander = st.sidebar.expander("Upload Variant List", expanded=False)

#         # with snp_list_upload_expander:

#         #     input_snp_file = st.file_uploader('  ')
        
#         # if input_snp_file:
#         #     input_snp_df = pd.read_csv(input_snp_file, header=None, names=['snp'])
#         #     input_snp_list = list(input_snp_df.loc[:,'snp'])
#         # ##### Add try except to check file format
#         # #### add check to see if anything in the input list is not in snp data
#         # else:
#         #     input_snp_list = None

#         if input_snp_list:
#             selected_snp = st.selectbox(
#             "Search by keywords for snpid, chr:pos, or chr:pos1-pos2", 
#             options= input_snp_list,
#             format_func=lambda x: 'Variant' if x == '' else x)
#         else:
#             # append empty string to front of snp list. if empty string chosen "select an option"
#             selected_snp = st.selectbox(
#                 "Search by keywords for snpid, chr:pos, or chr:pos1-pos2", 
#                 options=flagged_snps_list,
#                 format_func=lambda x: 'Variant' if x == '' else x)

#         if selected_snp:

#             geno_col = selected_snp + ".GType"
#             theta_col = selected_snp + ".Theta"
#             r_col = selected_snp + ".R"

#             gtypes = list(cb_df[geno_col].unique())

#             snp_for_plot_df = cb_df[['IID', theta_col, r_col, geno_col]].copy()
#             snp_for_plot_df.columns = snp_for_plot_df.columns.str.replace(f'{selected_snp}.', '', regex=False)

#             # df = snp_for_plot_df.copy()
#             # x_col, y_col = 'Theta', 'R'
#             # gtypes_list = list(df.GType.unique())

#             df = snp_for_plot_df.copy()
#             x_col, y_col, gtype_col = 'Theta', 'R', 'GType'
#             gtypes_list = (df.GType.unique())


#             gtypes_for_gmm = list(cb_df[geno_col].unique())
#             if 'NC' in gtypes_for_gmm:
#                 gtypes_for_gmm.remove('NC')
#             n_components = len(gtypes_for_gmm)

#             snp_for_gmm_df_temp = cb_df[['IID', theta_col, r_col, geno_col]].copy()
#             snp_for_gmm_df_temp.columns = snp_for_gmm_df_temp.columns.str.replace(f'{selected_snp}.','', regex=False)
#             snp_for_gmm_df = snp_for_gmm_df_temp.loc[~snp_for_gmm_df_temp['Theta'].isna() & ~snp_for_gmm_df_temp['R'].isna() & ~snp_for_gmm_df_temp['GType'].isna()].copy()

#             gmm_out = gtype_gmm(snp_theta_r_df=snp_for_gmm_df, n_components=n_components)

#             gmm_plot_df = gmm_out['X']
#             gmm_plot_df.loc[:,'GType'] = gmm_out['y_pred']
#             gmm_plot_df.loc[:,'IID'] = gmm_out['IID']
#             gmm = gmm_out['gmm']

#             cluster_means = {gtype:gmm_plot_df.loc[gmm_plot_df.GType==gtype,'Theta'].mean() for gtype in gmm_plot_df.GType.unique()}

#             gtype_order = ['AA','AB','BB']
#             gtypes_for_gmm_ordered = [gtype for gtype in gtype_order if gtype in gtypes_for_gmm]
#             ordered_cluster_means = dict(sorted(cluster_means.items(), key=lambda item: item[1]))
#             gmm_gtype_map = {list(ordered_cluster_means.keys())[i]:gtype for i, gtype in enumerate(gtypes_for_gmm_ordered)}

#             gmm_plot_df.loc[:,'gtype_out'] = gmm_plot_df.loc[:,'GType'].replace(gmm_gtype_map)
#             gmm_plot_df.loc[:,'original_gtype'] = snp_for_gmm_df.loc[:,'GType']
#             gmm_plot_df.loc[:,'snpid'] = selected_snp

#             reclustered_df = gmm_plot_df[['IID','snpid','gtype_out']].copy()
#             if "reclustered" not in st.session_state:
#                 st.session_state['reclustered'] = reclustered_df
#             else:
#                 st.session_state['reclustered'] = st.session_state['reclustered'].append(reclustered_df)

#             cluster_plot = plot_clusters(df, x_col, y_col, gtype_col=gtype_col, snpid=selected_snp)
#             cluster_fig = cluster_plot['fig']
#             xlim = cluster_plot['xlim']
#             ylim = cluster_plot['ylim']
#             recluster_fig = plot_hist_contour(df=gmm_plot_df, x_col='Theta', y_col='R', gtype_col='gtype_out', xlim=xlim, ylim=ylim)

#             left_column, right_column = st.columns([.8,1])
#             # left_column, right_column = st.columns([1.25,1])
            
#             with left_column:
#                 st.header("Default Clusters")
#                 left_exp = st.expander("Description", expanded=False)
#                 with left_exp:
#                     st.write('Genotypes assigned by the Illumina GenomeStudio Algorithm')

#                 # some padding to fix location of first plot
#                 # st.title('')
#                 # st.title('')
#                 st.write('')
#                 st.write('')
#                 st.write('')
#                 st.write('')
#                 st.write('')
#                 st.write('')
#                 st.write('')
#                 st.write('')
#                 # st.write('')
                
#                 st.plotly_chart(cluster_fig)

#             with right_column:
        
#                 st.header("DTi Clusters")
#                 right_exp = st.expander("Description", expanded=False)
#                 with right_exp:
#                     st.write('Genotypes assigned by DTi Enhanced Clustering Algorithm')

#                 st.plotly_chart(recluster_fig)
                
#             table_button = st.button('Show Table')
#             if table_button:
#                 st.table(snps_df.loc[snps_df.snpid == selected_snp].reset_index().drop(columns=['index']))

                
            
#             # fig = make_subplots(rows=1, cols=2)
#             # fig.add_trace(cluster_fig, row=1, col=1)
#             # fig.add_trace(recluster_fig, row=1, col=2)
#             # st.plotly_chart(fig)
#     st.sidebar.write("##")
#     export_expander = st.sidebar.expander("Download Reports", expanded=False)
        
#     with export_expander:

#         if 'reclustered' in st.session_state:
#             reclustered_df = st.session_state['reclustered']
#             if reclustered_df.shape[0] > 0:
#                 reclustered_csv = csv_convert_df(reclustered_df)
#                 st.download_button(
#                     label="Reclustered Variants",
#                     data=reclustered_csv,
#                     file_name='reclustered.csv',
#                     mime='text/csv'
#                 )

        
#         if 'flagged' in st.session_state:
#             flagged_out = st.session_state['flagged']

#             flagged_csv = csv_convert_df(flagged_out)
#             st.download_button(
#                 label="Flagged Variants",
#                 data=flagged_csv,
#                 file_name='flagged.csv',
#                 mime='text/csv'
#             )

#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')
#     st.sidebar.text(' ')

#     st.sidebar.image('img/DTI_logo_white_square-removebg-preview.png')