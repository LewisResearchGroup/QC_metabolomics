#!/usr/bin/env python
# coding: utf-8

# In[2]:

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


st.write('### Upload the results file generated on Mint')

results_file = st.file_uploader('results file')

try:
    s_st = SessionState.get(results = pd.read_csv(results_file))
    st.write('#### Your results file:')
    st.write(s_st.results.head())  
    
    
    st.write('#### indicate the intensity measurement')
    s_st.value_column = st.selectbox('select the intensity measurement \n', list(np.unique(s_st.results.columns)))
    
    
    st.write('#### indicate the growth control sample')
    s_st.gr_flag = st.text_input("growth control sample flag", 'ATCC')
    
    
    st.write('#### indicate the media control sample')
    s_st.me_flag = st.text_input("media control sample flag", 'MHPool')
    
    st.write('#### indicate the compound used for growing measurement')
    s_st.cp = st.selectbox('select the growing measurement compound \n', list(np.unique(s_st.results.peak_label)))
    
    st.write('#### indicate if the compound used for growing measurement is a consuming or a secreting one')
    s_st.flux = st.selectbox('select the compound interchange type \n', ['influx', 'eflux'])    
    
    st.write('#### set a threshold for the ratio between the growth control samples and the media control samples')
    s_st.threshold_0 = float(st.text_input("set a threshold for GControl/MSamples", '100'))
    
    div = np.mean(s_st.results[s_st.value_column][(s_st.results.peak_label == s_st.cp) & (s_st.results.ms_file.str.contains(s_st.gr_flag))]) / \
       np.mean(s_st.results[s_st.value_column][(s_st.results.peak_label == s_st.cp) & (s_st.results.ms_file.str.contains(s_st.me_flag))])
    
    
    st.write('the fraction of the compound between the growth media and the control control samples is:')
    st.write(1/div)
    
    if (s_st.flux == 'influx') & (div > 1/s_st.threshold_0):
        st.write('# Hey Sr, your controls didnt grow that well 😭😭😭😭😭')
        
    elif (s_st.flux == 'eflux') & (div < 1/s_st.threshold_0):
        st.write('# Hey Sr, your controls didnt grow that well 😭😭😭😭😭')    
    else:
        st.write('# Hey Sr, it looks like your controls grew well 🥳🥳🥳🥳🥳')
    
except:
    st.write('some point in your settings failed')



st.write(':heavy_minus_sign:' * 35)

st.write(':heavy_minus_sign:' * 35)
st.markdown('''# mz drift control section''')

def display_button2():
    
    display_instructions = st.selectbox('''Click here for explanation''' , ('Close', 'Show Explanation'))
    if display_instructions == 'Show Explanation':
        st.markdown('''in this section the column 'peak_mass_diff_50pc' is used to assess 
        the m/z drift a threshold value for the m/z drift should be provided. 
        This threshold indicates the m/z window size at which the half of the intensity is observed ''')
    
        
display_button2()

st.write('### Upload the results file generated on Mint, if they are the same that those in the previous section please continue ....')

results_file = st.file_uploader('results file for this section')

def get_percent_above(lista, thr):
    return 100*len(lista[abs(lista) > thr])/len(lista)




try:
    s_st.results = pd.read_csv(results_file)
except:
    pass


try:
    s_st = SessionState.get(results = pd.read_csv(results_file))
except:
    pass

try:
    st.write('#### Your results file:')
    st.write(s_st.results.head()) 

    st.write('#### set a threshold for the m/z drift')
    s_st.threshold_1 = float(st.text_input("set a threshold for m/z drift", '5'))
    
    st.write('searching ....')
    k = 0
    for metab in np.unique(s_st.results.peak_label):
        percent = get_percent_above(np.array(s_st.results.peak_mass_diff_50pc[s_st.results.peak_label == metab]), s_st.threshold_1)
#         print(percent)
        if percent > 10:
            k += 1
            if k == 1:
                st.write('### for the following compounds, more than 10 % of the samples have the m/z drift above the threshold: ')
            st.write(metab + ' :' + str(round(percent, 1)) + '%')

    
#     st.write('#### indicate the intensity measurement')
#     s_st.value_column = st.selectbox('select the intensity measurement \n', list(np.unique(s_st.results.columns)))
    
    
#     st.write('#### indicate the growth control sample')
#     s_st.gr_flag = st.text_input("growth control sample flag", 'ATCC')
    
    
#     st.write('#### indicate the media control sample')
#     s_st.me_flag = st.text_input("media control sample flag", 'MHPool')
    
#     st.write('#### indicate the compound used for growing measurement')
#     s_st.cp = st.selectbox('select the growing measurement compound \n', list(np.unique(s_st.results.peak_label)))
    
#     st.write('#### indicate if the compound used for growing measurement is a consuming or a secreting one')
#     s_st.flux = st.selectbox('select the compound interchange type \n', ['influx', 'eflux'])    
    
#     st.write('#### set a threshold for the ratio between the growth control samples and the media control samples')
#     s_st.threshold_0 = float(st.text_input("set a threshold for MSamples/GControl", '100'))
    
#     div = np.mean(s_st.results[s_st.value_column][(s_st.results.peak_label == s_st.cp) & (s_st.results.ms_file.str.contains(s_st.gr_flag))]) / \
#        np.mean(s_st.results[s_st.value_column][(s_st.results.peak_label == s_st.cp) & (s_st.results.ms_file.str.contains(s_st.me_flag))])
    
    
#     st.write('the fraction of the compound between the growth control samples and the control media is:')
#     st.write(div)
    
#     if (s_st.flux == 'influx') & (div > 1/s_st.threshold_0):
#         st.write('# Hey Sr, your controls didnt grow that well 😭😭😭😭😭')
        
#     elif (s_st.flux == 'eflux') & (div < 1/s_st.threshold_0):
#         st.write('# Hey Sr, your controls didnt grow that well 😭😭😭😭😭')    
#     else:
#         st.write('# Hey Sr, it looks like your controls grew well 🥳🥳🥳🥳🥳')
    
except:
    st.write('some point in your settings failed')


st.write(':heavy_minus_sign:' * 35)