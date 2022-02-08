import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from functools import reduce

from clusterbuster import csv_convert_df
from db_manhattan.manhattan import ManhattanPlot



def run():

    st.sidebar.title('DTi Genotype Analysis')
    reports_exp = st.sidebar.expander("Uploads", expanded=False)
    
    with reports_exp:
        report = st.file_uploader('IDATs')

        meta_data = st.file_uploader('Metadata')
        if report:
            df_qc = pd.read_csv('data/df_qc.csv')
            df_ancestry_counts = pd.read_csv('data/df_ancestry_counts.csv')
            df_ancestry_labels = pd.read_csv('data/df_ancestry_labels.csv')
            df_projected_pcs = pd.read_csv('data/df_projected_pcs.csv')
            df_ref_pcs = pd.read_csv('data/df_ref_pcs.csv')

            # df_ancestry_counts = pd.read_hdf(report, key='ancestry_counts')
            # df_ancestry_labels = pd.read_hdf(report, key='ancestry_labels')
            # df_confusion_matrix = pd.read_hdf(report, key='confusion_matrix')
            # df_new_samples_umap = pd.read_hdf(report, key='new_samples_umap')
            # df_projected_pcs = pd.read_hdf(report, key='projected_pcs')
            # df_ref_pcs = pd.read_hdf(report, key='ref_pcs')
            # df_ref_umap = pd.read_hdf(report, key='ref_umap')
            # df_total_umap = pd.read_hdf(report, key='total_umap')

            gwas = pd.read_csv('https://raw.githubusercontent.com/plotly/dash-bio-docs-files/master/manhattan_data.csv')

            st.session_state['df_qc'] = df_qc
            st.session_state['df_ancestry_counts'] = df_ancestry_counts
            st.session_state['df_ancestry_labels'] = df_ancestry_labels
            st.session_state['df_projected_pcs'] = df_projected_pcs
            st.session_state['df_ref_pcs'] = df_ref_pcs
            st.session_state['gwas'] = gwas
            st.session_state['pre_sample_n'] = 8000
            st.session_state['remaining_n'] = 8000

            qc_csv = csv_convert_df(df_qc)
            st.sidebar.download_button(
                label="Download QC report",
                data=qc_csv,
                file_name='qc.csv',
                mime='text/csv'
                )


    if 'df_qc' in st.session_state:
        df_qc = st.session_state['df_qc']
            
        # pre process
        # all sample prune
        df_2 = df_qc.query("level == 'sample'")
        df_2['sum'] = df_2.groupby('step')['pruned_count'].transform('sum')
        df_2 = df_2[['step','sum']]
        df_2 = df_2.drop_duplicates(subset=['step', 'sum'], keep='first')
        remaining_list = []

        start_n = st.session_state['remaining_n']
        for n in df_2['sum']:

            st.session_state['remaining_n']=st.session_state['remaining_n']-n
            rem = st.session_state['remaining_n']
            remaining_list.append(rem)


        remain_samples_x = [start_n] + remaining_list
        steps_y = ['pre_QC'] + list(df_2.step)
        funnel_df = pd.DataFrame({'remaining_samples':remain_samples_x,'step':steps_y})
        steps_dict = {
            'pre_QC': 'Pre-QC',
            'callrate_prune':'Call Rate Prune',
            'sex_prune': 'Sex Prune',
            'het_prune': 'Heterozygosity Prune',
            'related_prune': 'Relatedness Prune'
            }
        # funnel_df.loc[:,'step_name'] = funnel_df.replace({"step": steps_dict})
        funnel_df.loc[:,'step_name'] = funnel_df.loc[:,'step'].map(steps_dict)

        df_3 = df_qc.query("step == 'related_prune'")
        df_3 = df_3[['ancestry', 'pruned_count', 'metric']]

        df_3_related = df_3.query("metric == 'related_count'").reset_index(drop=True)
        df_3_related = df_3_related.rename(columns={'pruned_count': 'related_count'})
        df_3_related = df_3_related.drop('metric', 1)

        df_3_duplicated = df_3.query("metric == 'duplicated_count'").reset_index(drop=True)
        df_3_duplicated = df_3_duplicated.rename(columns={'pruned_count': 'duplicated_count'})
        df_3_duplicated = df_3_duplicated.drop('metric', 1)

        df_4 = pd.merge(df_3_related, df_3_duplicated, on="ancestry")
        ancestry_dict = {
            'FIN': 'Finnish (FIN)',
            'EAS': 'East Asian (EAS)',
            'AAC': 'African Admixted/Carribbean (AAC)',
            'AJ': 'Ashkenazi (AJ)',
            'SAS': 'South Asian (SAS)',
            'AMR': 'Native American/Latino (AMR)',
            'EUR': 'European (EUR)'
        }

        df_4.loc[:,'label'] = df_4.loc[:,'ancestry'].map(ancestry_dict)
        df_4.set_index('ancestry', inplace=True)

       
        #variant prune
        df_5 = df_qc.query("step == 'variant_prune'")
        df_5 = df_5[['ancestry', 'pruned_count', 'metric']]

        df_5_geno = df_5.query("metric == 'geno_removed_count'").reset_index(drop=True)
        df_5_geno = df_5_geno.rename(columns={'pruned_count': 'geno_removed_count'})
        df_5_geno = df_5_geno.drop('metric', 1)

        df_5_mis = df_5.query("metric == 'mis_removed_count'").reset_index(drop=True)
        df_5_mis = df_5_mis.rename(columns={'pruned_count': 'mis_removed_count'})
        df_5_mis = df_5_mis.drop('metric', 1)

        df_5_haplo = df_5.query("metric == 'haplotype_removed_count'").reset_index(drop=True)
        df_5_haplo = df_5_haplo.rename(columns={'pruned_count': 'haplotype_removed_count'})
        df_5_haplo = df_5_haplo.drop('metric', 1)

        df_5_hwe = df_5.query("metric == 'hwe_removed_count'").reset_index(drop=True)
        df_5_hwe = df_5_hwe.rename(columns={'pruned_count': 'hwe_removed_count'})
        df_5_hwe = df_5_hwe.drop('metric', 1)

        df_5_total = df_5.query("metric == 'total_removed_count'").reset_index(drop=True)
        df_5_total = df_5_total.rename(columns={'pruned_count': 'total_removed_count'})
        df_5_total = df_5_total.drop('metric', 1)

        data = [df_5_geno, df_5_mis, df_5_haplo, df_5_hwe, df_5_total]
        df_merged = reduce(lambda left,right: pd.merge(left,right,on=['ancestry'], how='outer'), data)
        df_merged.set_index('ancestry', inplace=True)

        df_6 = df_qc.loc[df_qc['pass'] == False]
        df_6 = df_6.reset_index(drop=True)


        # plotting
        #simple bar chart for callrate_prune and sex_prune
        funnel = px.funnel(
            funnel_df, 
            x='remaining_samples', 
            y='step_name', 
            height=600, width=550
            )

        #customize figure
        funnel.update_traces(
            marker_line_width=1.5, 
            opacity=0.8
            )
        funnel.update_layout(
        margin=dict(l=0, r=0, t=60, b=80),
        yaxis_title="QC Step"
    )
            # marker_color='rgb(158,202,225)', 
            # marker_line_color='rgb(8,48,107)',
        bar_3 = go.Figure(
            data=[
                go.Bar(y=df_4.label, x=df_4['related_count'], orientation='h', name="Related", base=0),
                go.Bar(y=df_4.label, x=-df_4['duplicated_count'], orientation='h', name="Duplicated", base=0)
                ])

        bar_3.update_layout(barmode='stack')

        bar_3.update_yaxes(
            ticktext=df_4.label,
            tickvals=df_4.label
        )

        bar_3.update_layout(
            autosize=False,
            height=600, width=700
        )

        bar_3.update_layout(
        margin=dict(l=0, r=0, t=60, b=80),
    )

        bar_6 = go.Figure(go.Bar(x=df_merged.index, y=df_merged['geno_removed_count'], name='Geno Removed Count'))
        bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['mis_removed_count'], name='Mis Removed Count'))
        bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['haplotype_removed_count'], name='Haplotype Removed Count'))
        bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['hwe_removed_count'], name='Hwe removed Count'))

        bar_6.update_layout(
            xaxis=dict(
                categoryorder='total descending',
                title='Ancestry',
                tickfont_size=14,
                
            ),
            yaxis=dict(
                title='Count',
                titlefont_size=16,
                tickfont_size=14,
            ),

            barmode='stack', 
            width=1200, height=600
        )


        if 'df_ancestry_counts' in st.session_state:
            df_ancestry_counts = st.session_state['df_ancestry_counts']
    
            #pie_chart for ancestry_counts
            df_ancestry_counts_percent = df_ancestry_counts.copy()
            df_ancestry_counts_percent['percent'] = 100 * df_ancestry_counts_percent['count']  / df_ancestry_counts_percent['count'].sum()
            pie_chart = px.pie(df_ancestry_counts_percent, values=df_ancestry_counts_percent.percent, names = df_ancestry_counts_percent.label)
            pie_chart.update_layout(showlegend=True,
            width=500,height=500
            )

        if 'df_projected_pcs' in st.session_state:
            df_projected_pcs = st.session_state['df_projected_pcs']
            # if selected_metrics_1 == 'Projected PCA':
            pcs13 = df_projected_pcs[['PC1', 'PC2', 'PC3']]
            
            pca_scatter = px.scatter_3d(
                pcs13, 
                x='PC1', 
                y='PC2', 
                z='PC3',
                color=df_projected_pcs['label'],
                width=600, height=600)
            pca_scatter.update_traces(marker_size = 4)


        # create app
        if 'remaining_n' in st.session_state:

            st.title('**GenoTools Quality Control**')
            
            st.header('QC Step 1: Sample Filtering')
            sample_exp = st.expander("Description", expanded=False)
            with sample_exp:
                st.write('\
                        1. Call Rate: Missingness per-individual > 0.02  \n\
                        2. Sex: F-coefficient < 0.25 is Female, F-coefficient > 0.75 is Male. 0.25 >= F <= 0.75 are outliers  \n\
                        3. Heterozygosity: Keep samples with -0.25 < F-coefficient < 0.25  \n\
                        4. Relatedness: Flag Samples with relatedness > 0.125 as cousins or more closely-related, Prune relatedness > 0.95 as Duplicates\
                ')

            left_col1, right_col1 = st.columns([1,1])
            with left_col1:
                st.header("**All Sample Filtering Counts**")
                st.plotly_chart(funnel)

            with right_col1:
                st.header("**Relatedness per Ancestry**")
                st.plotly_chart(bar_3)
            st.markdown('---')

            st.header('QC Step 2: Ancestry Estimation')
            anc_exp = st.expander("Description", expanded=False)
            with anc_exp:

                st.write('\
                    Labels assigned by the DTi ancestry algorithm trained on PC and UMAP-transformed 1000 Genomes genotypes\
                ')

            left_col2, right_col2 = st.columns([1.25,1])

            with right_col2:
                
                st.header("**Ancestry Distribution**")
                st.write(pie_chart)
                
            with left_col2:
                st.header("**Projected PCA**")
                st.plotly_chart(pca_scatter)
            st.markdown('---')

            st.header('QC Step 3: Variant Filtering')
            var_exp = st.expander("Description", expanded=False)
            with var_exp:
                st.write('\
                        1. Genotype Missingness <= 0.05  \n\
                        2. Case/Control Missingness: P > 1e-4  \n\
                        3. Haplotype Missingness: P > 1e-4  \n\
                        4. Hardy-Weinberg Equilibrium: P > 1e-4  \n\
                ')

            st.header("**Variant Filtering per Ancestry**")
            st.plotly_chart(bar_6)
            st.markdown('---')


            if df_6.shape[0]>0:
                st.markdown("**Failed Prune Steps**")
                st.table(df_6)
                st.markdown('---')
            
            

        st.text('  ')
        st.header('QC Step 4: Preliminary QC GWAS')
        if 'gwas' in st.session_state:
            gwas = st.session_state['gwas']

            gwas_exp = st.expander("Description", expanded=False)
            with gwas_exp:
                st.write("Summary statistics generated for largest ancestry group post-QC")
            
            st.header("**Manhattan Plot**")
            manhattan_fig = ManhattanPlot(
                dataframe = gwas
                )

            st.plotly_chart(manhattan_fig, use_container_width=True)

            st.header('**QQ-Plot**')
            
            st.image('img/qq_plot_final2.png')

        
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')
    st.sidebar.text(' ')

    st.sidebar.image('img/DTI_logo_white_square-removebg-preview.png')