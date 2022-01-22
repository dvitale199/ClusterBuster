from multiprocessing.connection import wait
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from sklearn.mixture import GaussianMixture
import time

from ClusterBuster.clusterbuster import parse_report, calculate_maf

st.write("""
# ClusterBuster
This webapp investigates cluster quality for genotyped variants from the **NeuroBoosterArray**
""")

def display_df_slice(df, max_rows=20, **st_dataframe_kwargs):
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
    n_rows = len(df)
    if n_rows <= max_rows:
        # As a special case, display small dataframe directly.
        st.write(df)
    else:
        # Slice the DataFrame to display less information.
        start_row = st.slider('Start row', 0, n_rows - max_rows)
        end_row = start_row + max_rows
        df = df[start_row:end_row]

        # Reindex Numpy arrays to make them more understadable.
        if type(df) == np.ndarray:
            df = pd.DataFrame(df)
            df.index = range(start_row,end_row)

        # Display everything.
        st.dataframe(df, **st_dataframe_kwargs)
        st.text('Displaying rows %i to %i of %i.' % (start_row, end_row - 1, n_rows))


# step 1: accept user input for report file

expander1 = st.sidebar.expander("Click Here to Upload a GenomeStudio SNP Report", expanded=False)
with expander1:
    st.subheader('Import GenomeStudio SNP Report')
    reportfile = st.file_uploader('')

    if reportfile:

        maf_threshold = st.slider('MAF Threshold for Flagging', value=0.010, step=0.005, min_value=0.000, max_value=1.0)
        gentrain_threshold = st.slider('MAF Threshold for Flagging', value=0.500, step=0.005, min_value=0.000, max_value=1.0)
        
        process_report_button = st.button("Process Report")

        if process_report_button:
            report_in = pd.read_csv(reportfile, engine='c', dtype={'Chr':str, 'position':int})
            report = parse_report(report_in, flag_maf=maf_threshold, flag_gentrain=gentrain_threshold)

            if report not in st.session_state:
                st.session_state['report'] = report
                st.session_state['update'] = False

            if report in st.session_state:
                st.session_state['report'] = report
                st.session_state['update'] = True

            with st.spinner('Calculating MAFs...'):
                time.sleep(2)
            with st.spinner('Flagging Variants for low MAF/GenCall...'):
                time.sleep(2)

            st.success('Done! Your file has been read and snps have been flagged for low MAF and GenTrain Score')

if 'report' in st.session_state.keys():
    
    report = st.session_state['report']
    cb_df = report['clusterbuster_df']
    snps_df = report['flagged_snps']

    # flagged dfs
    total_flagged = snps_df[snps_df.maf_flag | snps_df.gentrain_flag].reset_index()
    maf_flagged = snps_df.loc[snps_df.maf_flag].reset_index()
    gentrain_flagged = snps_df.loc[snps_df.gentrain_flag].reset_index()
    snps_list = snps_df.snpid.unique()
    flagged_list = total_flagged.snpid.unique()
    all_flagged = snps_df.loc[snps_df.snpid.isin(flagged_list)]

    expander2 = st.sidebar.expander("View Flagged SNPs", expanded=False)
    with expander2:
        
        maf_flagged_checkbox = st.checkbox('MAF', value=False)
        gentrain_flagged_checkbox = st.checkbox('GenTrain', value=False)

    if maf_flagged_checkbox:
        st.header('Flagged SNPs')
        # st.table(maf_flagged)
        display_df_slice(maf_flagged)
    if gentrain_flagged_checkbox:
        st.header('Flagged SNPs')
        display_df_slice(gentrain_flagged)
    if maf_flagged_checkbox and gentrain_flagged_checkbox:
        st.header('Flagged SNPs')
        display_df_slice(all_flagged)


    # st.header('Flagged SNPs')
    # st.table(snps_df.loc[snps_df.snpid.isin(flagged_list)])

        # selected_snp = st.selectbox("Select Flagged SNP to Plot", options=flagged_list)




    # expander2 = st.sidebar.expander("Search SNP", expanded=False)
    # with expander2:
    #     selected_snp = st.selectbox("Search by keywords for snpid, chr:pos, or chr:pos1-pos2", options=snps_list)
    # st.header('Selected SNP')
    # st.table(snps_df.loc[snps_df.snpid == selected_snp])







    
# step 3: search snps
# 4 different methods: 
#   1. by snp id
#   2. by chr:pos:a1:a2
#   3. by region --> chr:pos-pos
#   4. upload list of snps (1 column of snpids)
# for now, we will use this test list
# test_snps = pd.read_csv('data/testing_snps_of_interest.csv')
# SNP_list = test_snps.snp
# if report:


# st.sidebar.header('User Input Parameters, please type the name of your SNP of interest')





# step 4: return searched snps in table with id, chr, pos, gentrain, maf, frac{a,t,g,c} which can be selected to show plot
# add a highlight and click function to table. when clicked, return snpid to generate plot


# step 5: generate interactive clusterplot
#   1. clusters
#   2. labels for genotype
#   3. individual info --> ID, Theta, R, genotype
#   4. legend


# step 6: recluster using gmm
# same as step 5 plus ellipses + updated genotypes

# step 7: export cluster

# step 8: export updated snp report











# test_report = read_report('data/chr22_report.csv', flag_maf=0.01, flag_gencall=0.5)
# # SNP_list = test_report['flagged_snps'].snpid
# to_plot_df = test_report['clusterbuster_df']
# snps_df = test_report['flagged_snps']
# flagged = snps_df[snps_df.maf_flag & snps_df.gencall_flag]
# # SNP_list = flagged.snpid
# test_snps = pd.read_csv('data/testing_snps_of_interest.csv')
# SNP_list = test_snps.snp





# """# Cluster Plot."""





# """# Re-Cluster Plot."""
# def gmm_recluster(plot_df):

#     gtypes_for_gmm = list(plot_df[geno_col].unique())
#     gtypes_for_gmm.remove('NC')
#     n_components = len(gtypes_for_gmm)

#     X = plot_df[[theta_col,r_col]].copy()
#     gmm = GaussianMixture(n_components=n_components, covariance_type= "diag", random_state = 10).fit(X)
#     labels = gmm.predict(X)
#     X.loc[:,'predicted_label']=labels
#     # plt.figure(figsize=(9,7))
#     recluster_plot = sns.relplot(data=X, 
#                     x=theta_col,
#                     y=r_col, 
#                     hue="predicted_label")
#     recluster_plot.set(xlim=(-0.1, 1.1))
#     recluster_plot.set(ylim=(-0.1, 3.1))
#     st.pyplot(recluster_plot)

#     # X = plot_df[[theta_col,r_col]].copy()
#     # gmm = GaussianMixture(n_components=3, covariance_type= "full", random_state = 123).fit(X)
#     # labels = gmm.predict(X)
#     # X.loc[:,'predicted_label']=labels

#     # recluster_plot = sns.relplot(data=X, x=theta_col, y=r_col, hue="predicted_label")
#     # recluster_plot.set(xlim=(-0.1, 1.1))
#     # recluster_plot.set(ylim=(-0.1, 3.1))
#     # st.pyplot(recluster_plot)

# if st.button('RECLUSTER'):
#     gmm_recluster(to_plot_df)
#     # st.write('result: %s' % result)