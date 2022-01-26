import streamlit as st
# Imports here.
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import io

from clusterbuster import csv_convert_df, process_cnv_reports


def run():

    st.sidebar.title('DTi Genotype Analysis')
    bim_header = ["CHR","RS","CM","BP","A1","A2"]
    reports_exp = st.sidebar.expander("Upload CNV Reports", expanded=False)
    # st.session_state['plot_df'] = None
    # st.session_state['sample_id'] = None
    # st.session_state['gene_label'] = None

    with reports_exp:
        baf = st.file_uploader('BAF')
        lrr = st.file_uploader('LRR')
        bim = st.file_uploader('BIM')

        if baf and lrr and bim:
            BAF_temp = pd.read_csv(baf, engine='c', delim_whitespace=True)
            LRR_temp = pd.read_csv(lrr, engine='c', delim_whitespace=True)
            BIM = pd.read_csv(bim, engine='c', delim_whitespace=True, names=bim_header)
            st.session_state['BAF_temp'] = BAF_temp
            st.session_state['LRR_temp'] = LRR_temp
            st.session_state['BIM'] = BIM

    if all(df in st.session_state for df in ['BAF_temp','LRR_temp','BIM']):
        BAF_temp = st.session_state['BAF_temp']
        LRR_temp = st.session_state['LRR_temp']
        BIM = st.session_state['BIM']

        if BAF_temp.shape[0]>0 and LRR_temp.shape[0]>0 and BIM.shape[0]>0:
            gene_list = ['SNCA','SNCA','SNCA'] ####### REMOVE when we have full data#########
            sample_ids = [1029339,1056333,1398635] ####### REMOVE when we have full data #########
            cnv_type = ['Deletion','Duplication','Complex']
            # sample_ids = list(BAF_temp.IID.unique())
            genes = set(gene_list)
            sample_id = st.sidebar.selectbox('Sample ID', sample_ids)
            gene_label = st.sidebar.selectbox('Gene', genes) 
            process_report_button = st.sidebar.button("Analyze CNVs")

            if process_report_button:

                plot_df = process_cnv_reports(BAF_temp, LRR_temp, BIM, sample_id)

                # cnv_table = pd.DataFrame(
                #     {
                #         'sample_id':sample_ids,
                #         'gene':gene_list,
                #         'cnv_type':cnv_type
                #         }
                #     )
                # st.table(cnv_table) 
                ####### This data comes from reference later on
                chromosome = 4
                gene_start = 90645250
                gene_end = 90759466
                buffer = 1000000
                
                BAF_title = "Regional BAF Distribution"
                low_X = gene_start - buffer
                high_X = gene_end + buffer
                BAF_fig = px.scatter(plot_df, x='BP', y='BAF', color='LRR', title=BAF_title, color_continuous_scale='IceFire')
                BAF_fig.update_xaxes(range=[low_X, high_X])

                BAF_fig.add_shape(type="line",
                    x0=gene_start, y0=0.5, x1=gene_end, y1=0.5,
                    line=dict(color="Black",width=3)
                )

                annot_x = (gene_end + gene_start)/2
                annotation = {
                    # x -> location for x
                    'x': annot_x,
                    # y -> location for y
                    'y': 0.55,
                    'text': gene_label,  # text
                    'showarrow': True,  # would you want to see arrow
                    'arrowhead': 3,  # which type for arrowhead
                    'font': {'size': 10, 'color': 'black'}  # font style
                }

                BAF_fig.add_annotation(annotation)
                BAF_fig.update_layout()

                LRR_title = "Regional LRR Distribution"
                low_X = gene_start - buffer
                high_X = gene_end + buffer
                LRR_fig = px.scatter(plot_df, x='BP', y='LRR', color='BAF', title=LRR_title, color_continuous_scale='Twilight')
                LRR_fig.update_xaxes(range=[low_X, high_X])

                LRR_fig.add_shape(type="line",
                    x0=gene_start, y0=0.0, x1=gene_end, y1=0.0,
                    line=dict(color="Black",width=3)
                )

                annot_x = (gene_end + gene_start)/2
                annotation = {
                    # x -> location for x
                    'x': annot_x,
                    # y -> location for y
                    'y': 0.1,
                    'text': gene_label,  # text
                    'showarrow': True,  # would you want to see arrow
                    'arrowhead': 3,  # which type for arrowhead
                    'font': {'size': 10, 'color': 'black'}  # font style
                }

                LRR_fig.add_annotation(annotation)
                LRR_fig.update_layout()
                
                cnv_table = pd.DataFrame(
                    {'sample_id':sample_ids,
                    'gene':gene_list,
                    'cnv_type':cnv_type}
                    )
                cnv = cnv_table[(cnv_table.sample_id==sample_id) & (cnv_table.gene==gene_label)].reset_index().cnv_type[0]
                
                
                
                st.markdown("""
                <style>
                .big-font {
                    font-size:20px !important;
                }
                </style>
                """, unsafe_allow_html=True)

                st.markdown(f'<p class="big-font">{cnv} at {gene_label} on Chromosome {chromosome} +/- 1 MB</p>', unsafe_allow_html=True)
                
                # st.header(f'{cnv} at {gene_label} on Chromosome {chromosome} +/- 1 MB')
                st.plotly_chart(BAF_fig)
                st.plotly_chart(LRR_fig)
                

                download_exp = st.sidebar.expander("Download Reports", expanded=False)
                with download_exp:
                    

                    cnv_csv = csv_convert_df(cnv_table)
                    st.download_button(
                        label="CNVs",
                        data=cnv_csv,
                        file_name='cnvs.csv',
                        mime='text/csv'
                        )

                    buffer1 = io.StringIO()
                    BAF_fig.write_html(buffer1, include_plotlyjs='cdn')
                    html_bytes = buffer1.getvalue().encode()

                    st.download_button(
                        label='BAF Plot',
                        data=html_bytes,
                        file_name='BAF_fig.html',
                        mime='text/html'
                    )
                    
                    buffer2 = io.StringIO()
                    LRR_fig.write_html(buffer2, include_plotlyjs='cdn')
                    html_bytes = buffer2.getvalue().encode()

                    st.download_button(
                        label='LRR Plot',
                        data=html_bytes,
                        file_name='LRR_fig.html',
                        mime='text/html'
                    )
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
