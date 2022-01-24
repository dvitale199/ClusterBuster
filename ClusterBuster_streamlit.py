from multiprocessing.connection import wait
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from sklearn.mixture import GaussianMixture
import time
from PIL import Image

from clusterbuster import parse_report, view_table_slice, plot_clusters, gtype_gmm, plot_gmm, csv_convert_df
icon = Image.open('img/DTI_logo_white_square-removebg-preview.png')

st.set_page_config(
     page_title="DTi Genotype Analysis",
     page_icon=icon,
     layout="wide",
     initial_sidebar_state="expanded",
     menu_items={
         'About': "# This is a header. DTi Genotype Analysis: Investigate cluster quality for genotyped variants from the **NeuroBoosterArray**"
     }
 )

st.title('DTi Genotype Analysis')

# st.write("""
# # ClusterBuster
# This webapp investigates cluster quality for genotyped variants from the **NeuroBoosterArray**
# """)

# step 1: accept user input for report file
# expander1 = st.sidebar.expander("Click Here to Upload a GenomeStudio SNP Report", expanded=False)
# with expander1:
st.sidebar.subheader('Upload GenomeStudio SNP Report')
reportfile = st.sidebar.file_uploader('')

if reportfile:

    maf_threshold = st.sidebar.slider('MAF Threshold for Flagging', value=0.010, step=0.005, min_value=0.000, max_value=1.0)
    gentrain_threshold = st.sidebar.slider('GenTrain Threshold for Flagging', value=0.500, step=0.005, min_value=0.000, max_value=1.0)
    
    process_report_button = st.sidebar.button("Process Report")

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
    snps_list = list(snps_df.snpid.unique())

    # flagged dfs
    total_flagged = snps_df[snps_df.maf_flag | snps_df.gentrain_flag].reset_index()
    both_flagged = snps_df[snps_df.maf_flag & snps_df.gentrain_flag].reset_index()
    maf_flagged = snps_df.loc[snps_df.maf_flag].reset_index()
    gentrain_flagged = snps_df.loc[snps_df.gentrain_flag].reset_index()
    reclustered_df = pd.DataFrame()

    snps_list = list(snps_df.snpid.unique())

    flagged_list = list(total_flagged.snpid.unique())
    all_flagged = snps_df.loc[snps_df.snpid.isin(flagged_list)]

    expander2 = st.sidebar.expander("Filter Flagged SNPs", expanded=False)
    with expander2:
        
        maf_flagged_checkbox = st.checkbox('MAF', value=False)
        gentrain_flagged_checkbox = st.checkbox('GenTrain', value=False)
        all_flagged_checkbox = st.checkbox('All', value=False)
    
    if maf_flagged_checkbox and not gentrain_flagged_checkbox:
        st.header('Flagged SNPs')
        view_table_slice(maf_flagged)
        all_flagged_checkbox = False

    if gentrain_flagged_checkbox and not maf_flagged_checkbox:
        st.header('Flagged SNPs')
        view_table_slice(gentrain_flagged)
        all_flagged_checkbox = False

    if maf_flagged_checkbox and gentrain_flagged_checkbox:
        st.header('Flagged SNPs')
        view_table_slice(both_flagged)
        all_flagged_checkbox = False

    if all_flagged_checkbox:
        st.header('Flagged SNPs')
        view_table_slice(total_flagged)
        gentrain_flagged_checkbox = False
        maf_flagged_checkbox = False

    expander3 = st.sidebar.expander("Search SNPs", expanded=False)

    with expander3:

        st.write('Upload a list of Variants:')
        input_snp_file = st.file_uploader('Variant List Upload')
        if input_snp_file:
            input_snp_df = pd.read_csv(input_snp_file, header=None, names=['snp'])
            input_snp_list = list(input_snp_df.loc[:,'snp'])
        ##### Add try except to check file format
        #### add check to see if anything in the input list is not in snp data
        else:
            input_snp_list = None

        if input_snp_list:
            selected_snp = st.selectbox(
            "Search by keywords for snpid, chr:pos, or chr:pos1-pos2", 
            options=[''] + input_snp_list,
            format_func=lambda x: 'Select an option' if x == '' else x)
        else:
            # append empty string to front of snp list. if empty string chosen "select an option"
            selected_snp = st.selectbox(
                "Search by keywords for snpid, chr:pos, or chr:pos1-pos2", 
                options=[''] + snps_list,
                format_func=lambda x: 'Select an option' if x == '' else x)

    if selected_snp:

        st.table(snps_df.loc[snps_df.snpid == selected_snp].reset_index().drop(columns=['index']))

        geno_col = selected_snp + ".GType"
        theta_col = selected_snp + ".Theta"
        r_col = selected_snp + ".R"

        gtypes = list(cb_df[geno_col].unique())

        snp_for_plot_df = cb_df[['IID', theta_col, r_col, geno_col]].copy()
        snp_for_plot_df.columns = snp_for_plot_df.columns.str.replace(f'{selected_snp}.', '', regex=False)

        # df = snp_for_plot_df.copy()
        # x_col, y_col = 'Theta', 'R'
        # gtypes_list = list(df.GType.unique())

        df = snp_for_plot_df.copy()
        x_col, y_col, gtype_col = 'Theta', 'R', 'GType'
        gtypes_list = (df.GType.unique())

        st.plotly_chart(plot_clusters(df, x_col, y_col, gtype_col=gtype_col, snpid=selected_snp))

        gtypes_for_gmm = list(cb_df[geno_col].unique())
        gtypes_for_gmm.remove('NC')
        n_components = len(gtypes_for_gmm)

        snp_for_gmm_df_temp = cb_df[['IID', theta_col, r_col, geno_col]].copy()
        snp_for_gmm_df_temp.columns = snp_for_gmm_df_temp.columns.str.replace(f'{selected_snp}.','', regex=False)
        snp_for_gmm_df = snp_for_gmm_df_temp.loc[~snp_for_gmm_df_temp['Theta'].isna() & ~snp_for_gmm_df_temp['R'].isna() & ~snp_for_gmm_df_temp['GType'].isna()].copy()
        

        recluster = st.button('Recluster Selected Variant')
        

        if recluster:

            gmm_out = gtype_gmm(snp_theta_r_df=snp_for_gmm_df, n_components=n_components)

            gmm_plot_df = gmm_out['X']
            gmm_plot_df.loc[:,'GType'] = gmm_out['y_pred']
            gmm_plot_df.loc[:,'IID'] = gmm_out['IID']
            gmm = gmm_out['gmm']

            cluster_means = {gtype:gmm_plot_df.loc[gmm_plot_df.GType==gtype,'Theta'].mean() for gtype in gmm_plot_df.GType.unique()}

            gtype_order = ['AA','AB','BB']
            gtypes_for_gmm_ordered = [gtype for gtype in gtype_order if gtype in gtypes_for_gmm]
            ordered_cluster_means = dict(sorted(cluster_means.items(), key=lambda item: item[1]))
            gmm_gtype_map = {list(ordered_cluster_means.keys())[i]:gtype for i, gtype in enumerate(gtypes_for_gmm_ordered)}

            gmm_plot_df.loc[:,'gtype_out'] = gmm_plot_df.loc[:,'GType'].replace(gmm_gtype_map)
            gmm_plot_df.loc[:,'original_gtype'] = snp_for_gmm_df.loc[:,'GType']
            gmm_plot_df.loc[:,'snpid'] = selected_snp
            reclustered_df = reclustered_df.append(gmm_plot_df[['IID','snpid','gtype_out']])

            st.plotly_chart(plot_gmm(df=gmm_plot_df, x_col='Theta', y_col='R', gtype_col='gtype_out', gmm=gmm, snpid=selected_snp, n_std=5))

    else:
        st.warning('Please Select a SNP')
    

    if reclustered_df.shape[0] > 0:
        reclustered_csv = csv_convert_df(reclustered_df)
        st.sidebar.download_button(
            label="Download Reclustered Variants",
            data=reclustered_csv,
            file_name='reclustered.csv',
            mime='text/csv'
        )

        
    # missing_df = pd.DataFrame()

    # # Generate output report
    # for snp in snps_list:
    #     geno_col = snp + ".GType"
    #     theta_col = snp + ".Theta"
    #     r_col = snp + ".R"    

    #     gtypes = list(cb_df[geno_col].unique())

    #     snp_for_plot_df = cb_df[['IID', theta_col, r_col, geno_col]].copy()
    #     snp_for_plot_df.columns = snp_for_plot_df.columns.str.replace(f'{snp}.', '', regex=False)

    #     df = snp_for_plot_df.copy()
    #     x_col, y_col = 'Theta', 'R'
    #     gtypes_list = list(df.GType.unique())

    #     df = snp_for_plot_df.copy()
    #     x_col, y_col, gtype_col = 'Theta', 'R', 'GType'
    #     gtypes_list = (df.GType.unique())

    #     gtypes_for_gmm = list(cb_df[geno_col].unique())
    #     # if 'NC' in 
    #     # gtypes_for_gmm.remove('NC')
    #     n_components = len(gtypes_for_gmm)

    #     snp_for_gmm_df_temp = cb_df[['IID', theta_col, r_col, geno_col]].copy()
    #     snp_for_gmm_df_temp.columns = snp_for_gmm_df_temp.columns.str.replace(f'{snp}.','', regex=False)
    #     # print(snp_for_gmm_df_temp.head())
    #     snp_for_gmm_df = snp_for_gmm_df_temp.loc[~snp_for_gmm_df_temp['Theta'].isna() & ~snp_for_gmm_df_temp['R'].isna() & ~snp_for_gmm_df_temp['GType'].isna()].copy()

    #     # snp_for_gmm_df = snp_for_gmm_df_temp.loc[~snp_for_gmm_df_temp.Theta.isna() & ~snp_for_gmm_df_temp.R.isna() & ~snp_for_gmm_df_temp.GType.isna()].copy()
    #     # snp_for_gmm_df.columns = snp_for_gmm_df.columns.str.replace(f'{snp}.','', regex=False)

    #     # handle missing
        

    #     missing_snp = snp_for_gmm_df.loc[snp_for_gmm_df.Theta.isna() | snp_for_gmm_df.R.isna() | snp_for_gmm_df.GType.isna()].copy()
    #     if missing_snp.shape[0] > 0:
    #         missing_snp.loc[:,'snpid'] = snp
    #         missing_df = missing_df.append(missing_snp)

    # missing_csv = csv_convert_df(missing_df)
    # all_flagged_csv = csv_convert_df(all_flagged)



    # expander4 = st.sidebar.expander("Downloads", expanded=False)
    # with expander4:
    #     st.sidebar.download_button(
    #         label="Download samples with missing theta|r|genotype",
    #         data=missing_csv,
    #         file_name='missing.csv',
    #         mime='text/csv'
    #     )

    #     st.sidebar.download_button(
    #         label="Download Flagged Variants File",
    #         data=all_flagged_csv,
    #         file_name='flagged.csv',
    #         mime='text/csv'
    #     )



        

        










    
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