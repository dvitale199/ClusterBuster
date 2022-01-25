import streamlit as st
import pandas as pd
# import hdfStore
import tempfile
import io
import deepdish as dd
import plotly.express as px
import plotly.graph_objects as go
from functools import reduce
import dash_bio as dashbio

def plot_manhattan(df):
    fig = dashbio.ManhattanPlot(
        dataframe = df
    )
    # figure.update_layout(showlegend=True,
    #     width=1200,
    #     height=500)

    # st.write(fig)
    fig.write_html("img/manhattan.html")

def run():
    st.sidebar.title('DTi Genotype Analysis')
    reports_exp = st.sidebar.expander("Upload GenoTools Reports", expanded=False)

    with reports_exp:
        report = st.file_uploader('')
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

            
            selected_metrics = st.sidebar.selectbox(label="Prune selection", options=['All Sample Prune', 'Related Prune', 'Variant Prune'])
            selected_metrics_1 = st.sidebar.selectbox(label="PCA selection", options=['Reference PCA', 'Projected PCA'])
            
            st.session_state['selected_metrics'] = selected_metrics
            st.session_state['selected_metrics_1'] = selected_metrics_1
            st.session_state['df_qc'] = df_qc
            st.session_state['df_ancestry_counts'] = df_ancestry_counts
            st.session_state['df_ancestry_labels'] = df_ancestry_labels
            st.session_state['df_projected_pcs'] = df_projected_pcs
            st.session_state['df_ref_pcs'] = df_ref_pcs
            st.session_state['gwas'] = gwas

    if 'selected_metrics' in st.session_state:
        selected_metrics = st.session_state['selected_metrics']
        if 'df_qc' in st.session_state:
            df_qc = st.session_state['df_qc']
            if selected_metrics == 'All Sample Prune':
                #all sample prune
                df_2 = df_qc.query("level == 'sample'")
                df_2['sum'] = df_2.groupby('step')['pruned_count'].transform('sum')
                df_2 = df_2[['step','sum']]
                df_2 = df_2.drop_duplicates(subset=['step', 'sum'], keep='first')

                print(df_2.head())

                #simple bar chart for callrate_prune and sex_prune

                bar_2 = px.bar(df_2, x='step', y='sum', text='sum',
                            hover_data=['step','sum'], 
                            labels={'step':'Pruning Step', 'sum':'Count'}, height=500, width=800)

                #customize figure
                bar_2.update_traces(marker_color='rgb(158,202,225)', marker_line_color='rgb(8,48,107)',
                                marker_line_width=1.5, opacity=0.6)
                st.markdown("**All Sample Prune Count**")
                st.write(bar_2)


            if selected_metrics == 'Related Prune':
                #related prune
                df_3 = df_qc.query("step == 'related_prune'")
                df_3 = df_3[['ancestry', 'pruned_count', 'metric']]

                df_3_related = df_3.query("metric == 'related_count'").reset_index(drop=True)
                df_3_related = df_3_related.rename(columns={'pruned_count': 'related_count'})
                df_3_related = df_3_related.drop('metric', 1)

                df_3_duplicated = df_3.query("metric == 'duplicated_count'").reset_index(drop=True)
                df_3_duplicated = df_3_duplicated.rename(columns={'pruned_count': 'duplicated_count'})
                df_3_duplicated = df_3_duplicated.drop('metric', 1)

                df_4 = pd.merge(df_3_related, df_3_duplicated, on="ancestry")
                df_4.set_index('ancestry', inplace=True)

                bar_3 = go.Figure(data=[
                go.Bar(y=df_4.index, x=df_4['related_count'], orientation='h', name="Related Count", base=0),
                go.Bar(y=df_4.index, x=-df_4['duplicated_count'], orientation='h', name="Duplicated Count", base=0)
                ])

                bar_3.update_layout(
                barmode='stack')

                bar_3.update_yaxes(
                    ticktext=df_4.index,
                    tickvals=df_4.index
                )
                st.markdown("**Related Prune Count per Ancestry**")
                st.write(bar_3)

            if selected_metrics == 'Variant Prune':
                #vaient prune
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

                bar_6 = go.Figure(go.Bar(x=df_merged.index, y=df_merged['geno_removed_count'], name='Geno Removed Count'))
                bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['mis_removed_count'], name='Mis Removed Count'))
                bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['haplotype_removed_count'], name='Haplotype Removed Count'))
                bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['hwe_removed_count'], name='Hwe removed Count'))

                bar_6.update_layout(
                    xaxis=dict(
                        title='Ancestry',
                        tickfont_size=14,
                    ),
                    yaxis=dict(
                        title='Count',
                        titlefont_size=16,
                        tickfont_size=14,
                    )
                )

                bar_6.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'},width=800,height=600)
                st.markdown("**Variant Prune Count per Ancestry**")
                st.write(bar_6)

            df_6 = df_qc.loc[df_qc['pass'] == False]
            df_6 = df_6.reset_index(drop=True)
            if df_6.shape[0]>0:
                st.markdown("**Failed Prune Steps**")
                st.table(df_6)
            
            left_column, right_column = st.columns([1.25,1])
                
            with right_column:
                if 'df_ancestry_counts' in st.session_state:
                    df_ancestry_counts = st.session_state['df_ancestry_counts']
            
                    #pie_chart for ancestry_counts
                    df_ancestry_counts_percent = df_ancestry_counts.copy()
                    df_ancestry_counts_percent['percent'] = 100 * df_ancestry_counts_percent['count']  / df_ancestry_counts_percent['count'].sum()
                    pie_chart = px.pie(df_ancestry_counts_percent, values=df_ancestry_counts_percent.percent, names = df_ancestry_counts_percent.label)
                    pie_chart.update_layout(showlegend=True,
                        width=550,
                        height=550)
                    st.markdown("**Ancestry Distribution**")
                    st.write(pie_chart)

            with left_column:
                if 'selected_metrics_1' in st.session_state:
                    selected_metrics_1 = st.session_state['selected_metrics_1']
                    if 'df_ref_pcs' in st.session_state:
                        df_ref_pcs = st.session_state['df_ref_pcs']
                        if selected_metrics_1 == 'Reference PCA':
                            X = df_ref_pcs[['PC1', 'PC2', 'PC3']]

                            fig = px.scatter_3d(
                                X, 
                                x='PC1', 
                                y='PC2', 
                                z='PC3',
                                color=df_ref_pcs['label'])
                            
                            st.markdown("**Reference PCA**")
                            st.plotly_chart(fig)
                    
                    if 'df_projected_pcs' in st.session_state:
                        df_projected_pcs = st.session_state['df_projected_pcs']
                        if selected_metrics_1 == 'Projected PCA':
                            X = df_projected_pcs[['PC1', 'PC2', 'PC3']]
                            
                            fig = px.scatter_3d(
                                X, 
                                x='PC1', 
                                y='PC2', 
                                z='PC3',
                                color=df_projected_pcs['label'])
                            
                            st.markdown("**Projected PCA**")
                            st.plotly_chart(fig)

        #manhattan plot
    if 'gwas' in st.session_state:
        gwas = st.session_state['gwas']
        DATASET = gwas.groupby('CHR').apply(lambda u: u.head(50))
        DATASET = DATASET.droplevel('CHR').reset_index(drop=True)
        
        plot_manhattan(gwas)

    # threshold = st.slider('Threshold value', min_value = 1, max_value = 10, value =1)
    # state = st.session_state.get(position=0)
    # widget = st.empty()
    # ,on_change=plot_manhattan(df1, )