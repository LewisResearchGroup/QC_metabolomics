#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import os
import pandas as pd
import datetime
import numpy as np
import glob
import re
import streamlit as st
from stlibs import SessionState



# In[3]:


def download_link(object_to_download, download_filename, download_link_text):
    """
    Generates a link to download the given object_to_download.
    object_to_download (str, pd.DataFrame):  The object to be downloaded.
    download_filename (str): filename and extension of file. e.g. mydata.csv, some_txt_output.txt
    download_link_text (str): Text to display for download link.
    Examples:
    download_link(YOUR_DF, 'YOUR_DF.csv', 'Click here to download data!')
    download_link(YOUR_STRING, 'YOUR_STRING.txt', 'Click here to download your text!')
    """
    if isinstance(object_to_download,pd.DataFrame):
        object_to_download = object_to_download.to_csv(index=False)

    # some strings <-> bytes conversions necessary here
    b64 = base64.b64encode(object_to_download.encode()).decode()


# In[4]:


def display_button():
    display_instructions = st.selectbox('''Click here to see explanation''' , ('Close', 'Show Explanation'))
    if display_instructions == 'Show Explanation':
        st.markdown("""
         #### This app is meant to be used internally for running metabolomics quality control (QC)
         #### The QC pipeline is running in five phases:
         """)
        st.write('''         
         #####    1) Growth control
         #####    2) M/Z drift calculation
         #####    3) RT drift (not sure yet)    
         #####    4) Linear Ranges
         #####    5) Detection limits based on signal to noise ratio (LOD and LOQ)
         ''')
        
display_button()

st.write(':heavy_minus_sign:' * 35)
st.markdown('''# Growth control section''')


st.write('### Upload the results file and indicate the program used to generate it')

results_file = st.file_uploader('results file')



try:
    s_st = SessionState.get(results = pd.read_csv(results_file))
    st.write('#### Your results concentrations file:')
    st.write(s_st.results.head())  
    
    
    st.write('#### indicate the intensity measurement')
    s_st.value_column = st.selectbox('select the intensity measurement \n', list(np.unique(s_st.results.columns)))
    
    
    st.write('#### indicate the growth control sample')
    s_st.gr_flag = st.text_input("growth control sample flag", 'ATCC')
    
    
    st.write('#### indicate the media control sample')
    s_st.me_flag = st.text_input("media control sample flag", 'MHPool')
    
    st.write('#### indicate the consuming compound used for growing measurement')
    s_st.cp = st.selectbox('select the growing measurement compound \n', list(np.unique(s_st.results.peak_label)))
    
    
    st.write('#### set a threshold for the ratio between the growth control samples and the media control samples')
    s_st.threshold_0 = int(st.text_input("set a threshold for GControl/MSamples", '10'))
    
    div = np.mean(s_st.results[s_st.value_column][(s_st.results.peak_label == s_st.cp) & (s_st.results.ms_file.str.contains(s_st.gr_flag))]) / \
       np.mean(s_st.results[s_st.value_column][(s_st.results.peak_label == s_st.cp) & (s_st.results.ms_file.str.contains(s_st.me_flag))])
    
    st.write(div)
    if div > 1/s_st.threshold_0:
        st.write('# Hey Sr, your controls didnt grow that well')
    
except:
    st.write('some point in your settings failed')

# std_info = st.sidebar.file_uploader('Upload results file')

# st.write('enter the name of growth control samples flag')


st.write(':heavy_minus_sign:' * 35)

