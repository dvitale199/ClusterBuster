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

icon = Image.open('img/dti_gray_logo.png')

st.set_page_config(
     page_title="DTi Genotype Analysis",
     page_icon=icon,
     layout="wide",
     initial_sidebar_state="expanded",
     menu_items={
         'About': "# This is a header. DTi Genotype Analysis: Investigate cluster quality for genotyped variants from the **NeuroBoosterArray**"
     }
 )


st.markdown(
    '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@4.5.3/dist/css/bootstrap.min.css" integrity="sha384-TX8t27EcRE3e/ihU7zmQxVncDAy5uIKz4rEkgIXeMed4M0jlfIDPvg6uqKI2xXr2" crossorigin="anonymous">',
    unsafe_allow_html=True,)
query_params = st.experimental_get_query_params()
tabs = ["ClusterBuster", "PlusCNV", "GenoTools"]
if "tab" in query_params:
    active_tab = query_params["tab"][0]
else:
    active_tab = "ClusterBuster"

if active_tab not in tabs:
    st.experimental_set_query_params(tab="ClusterBuster")
    active_tab = "ClusterBuster"

li_items = "".join(
    f"""
    <li class="nav-item">
        <a class="nav-link{' active' if t==active_tab else ''}" href="/?tab={t}">{t}</a>
    </li>
    """
    for t in tabs
)
tabs_html = f"""
    <ul class="nav nav-tabs">
    {li_items}
    </ul>
"""

st.markdown(tabs_html, unsafe_allow_html=True)
st.markdown("<br>", unsafe_allow_html=True)

if active_tab == "ClusterBuster":
    cb.run()
elif active_tab == "PlusCNV":
    cnv.run()
elif active_tab == "GenoTools":
    gt.run()
else:
    st.error("Something has gone terribly wrong.")
