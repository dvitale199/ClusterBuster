import streamlit as st
# Imports here.
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os


def process_cnv_reports(BAF, LRR, BIM, sample_id):

    # BAF reduce and transpose.
    BAF_temp_reduced = BAF[BAF['V1'] == sample_id]
    BAF_temp_reduced.drop(columns=['V1','V2','V3','V4','V5','V6'], inplace=True)
    BAF_transposed = BAF_temp_reduced.transpose()
    BAF_transposed.columns = ["BAF"]

    # LRR reduce and transpose.
    LRR_temp_reduced = LRR[LRR['V1'] == sample_id]
    LRR_temp_reduced.drop(columns=['V1','V2','V3','V4','V5','V6'], inplace=True)
    LRR_transposed = LRR_temp_reduced.transpose()
    LRR_transposed.columns = ["LRR"]

    BAF_transposed.reset_index(drop=True, inplace=True)
    LRR_transposed.reset_index(drop=True, inplace=True)
    BIM.reset_index(drop=True, inplace=True)

    out_df = pd.concat([BAF_transposed, LRR_transposed, BIM], axis=1)

    return out_df

def run():
    
    st.title('PlusCNV')

    reports_exp = st.sidebar.expander("Upload CNV Reports", expanded=False)
    with reports_exp:
        baf = st.file_uploader('BAF')
        lrr = st.file_uploader('LRR')
        bim = st.file_uploader('BIM')

        process_report_button = st.sidebar.button("Analyze CNVs")

        if process_report_button:

    ##### TEST FILES ####
    # baf : "SNCA_gene_region_1_ukb_baf_chr4_v2_FINAL_reduced.txt"
    # lrr : "SNCA_gene_region_1_ukb_l2r_chr4_v2_FINAL_reduced.txt"
    # bim : "SNCA_gene_region1_ukb_bim_chr4_v2.txt"
    ###################################
    
            BAF_temp = pd.read_csv(baf, engine='c', delim_whitespace=True)
            LRR_temp = pd.read_csv(lrr, engine='c', delim_whitespace=True)
            bim_header = ["CHR","RS","CM","BP","A1","A2"]
            BIM = pd.read_csv(bim, engine='c', delim_whitespace=True, names=bim_header)

            plot_df = process_cnv_reports(BAF_temp, LRR_temp, BIM, sample_id)
            #### Check if files are read #####

            
            # BAF_temp = True ####### REMOVE #########
            # if BAF_temp:

            gene_list = ['SNCA'] ####### REMOVE #########
            sample_ids = [1029339] ####### REMOVE #########
            # sample_ids = list(BAF_temp.IID.unique())

            sample_id = st.sidebar.selectbox('Sample ID', sample_ids)
        
            gene_label = st.sidebar.selectbox('Gene', gene_list)
        


            # # Select the gene. This is where positions would be looked up. That can come later after demo.
            # # Also should include a window of X BP around gene as an option for viewing.
            chromosome = 4
            gene_start = 90645250
            gene_end = 90759466
            buffer = 1000000


            


        # def plot_BAF():
            BAF_title = "sample ID " + str(sample_id) + " @ " + gene_label + " on CHR " + str(chromosome) + " +/- " + str(buffer) + " BP."
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

            st.plotly_chart(BAF_fig)



            LRR_title = "sample ID " + str(sample_id) + " @ " + gene_label + " on CHR " + str(chromosome) + " +/- " + str(buffer) + " BP."
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

            st.plotly_chart(LRR_fig)