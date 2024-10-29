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
         #####    4) Peak intensity
         ''')
        
display_button()


st.write(':heavy_minus_sign:' * 35)
st.markdown('''# parameters drift control section''')



def display_button2():
    
    display_instructions = st.selectbox('''Click here for explanation''' , ('Close', 'Show Explanation'))
    if display_instructions == 'Show Explanation':
        st.markdown('''in this section you will be comparing the new generated data from a historical one to compare the retention time differences. 
        The analysis will be restricted to the standard samples of the same concentrations. The column 'peak_rt_of_max' of mint output will be used for such comparison.     Standard samples will be considered as those containing the Std flag ''')
    
        
display_button2()

st.write('### Upload the historical results file generated on Mint ....')
historical_file = st.file_uploader('historical results')


try:
    s_st = SessionState.get( historical_results = pd.read_csv(historical_file ))
except:
    pass

try:
    s_st.historical_results = pd.read_csv(historical_file)    
except:
    pass


try:
    st.write('#### Your historical file:')
    st.write(s_st.historical_results.head()) 
    
except:
    st.write('some point in your settings failed')

st.write('### Upload the newly generated results file ....')
results_file = st.file_uploader('current results')

try:
    s_st.results = pd.read_csv(results_file)
    s_st.results0 = s_st.results.copy()
except:
    pass



try:
    st.write('#### Your results file:')
    st.write(s_st.results.head()) 
    
except:
    st.write('some point in your settings failed')

try:
    lres = len(s_st.results)
    lhres = len(s_st.historical_results)
except:
    lres = 0
    lhres = 0
if (lres > 1) & (lhres > 1):
    
    ### selecting the standard samples only ###
    try: 
        st.write('#### indicate the Standard samples flag')
        s_st.std_flag = st.text_input("growth control sample flag", 'Std')
        s_st.historical_results = s_st.historical_results[s_st.historical_results.ms_file.str.contains(s_st.std_flag)]
        st.write('your historical data has ' + str(len(np.unique(s_st.historical_results.ms_file))) + ' ' + s_st.std_flag + ' samples')

    except:
        st.write('there is a problem with the standard sample selection in the historical results')
    
    try:
        s_st.results = s_st.results[s_st.results.ms_file.str.contains(s_st.std_flag)]
        st.write('your results data has ' + str(len(np.unique(s_st.results.ms_file))) + ' ' + s_st.std_flag + ' samples')
        
        if len(np.unique(s_st.results.ms_file)) == len(np.unique(s_st.historical_results.ms_file)):
            st.write('which matches the size of the historical results')
    except:
        st.write('there is a problem with the standard sample selection in the current results')
    
    ### numbering the standard samples by the file name ###
    try:
        s_st.historical_results['STDType'] = s_st.historical_results.ms_file.apply(lambda x: x.split('.')[0])
        s_st.historical_results.STDType = s_st.historical_results.STDType.apply(lambda x: int(x.split('Std')[-1]))
        
        s_st.results['STDType'] = s_st.results.ms_file.apply(lambda x: x.split('.')[0])
        s_st.results.STDType = s_st.results.STDType.apply(lambda x: int(x.split('Std')[-1]))
        
        st.write('there are ' + str(len(np.unique(s_st.historical_results.STDType))) + ' types of ' + s_st.std_flag + ' sample types in the historical results file')
        st.write('there are ' + str(len(np.unique(s_st.results.STDType))) + ' types of ' + s_st.std_flag + ' sample types in the results file')
        s_st.intersection_samples =  np.intersect1d(np.unique(s_st.historical_results.STDType), np.unique(s_st.results.STDType))
        st.write('with ' + str(len(s_st.intersection_samples)) + ' samples in the intersection')
    except:
        st.write('there is a problem with the standard sample numbering')
    
    ### getting compounds that exist in historical and current data ###
    try:
        st.write('there are ' + str(len(np.unique(s_st.historical_results.peak_label))) +' compounds in the historical data')
        st.write('there are ' + str(len(np.unique(s_st.results.peak_label))) +' compounds in the results data')
        s_st.intersection_compounds = np.intersect1d(np.unique(s_st.historical_results.peak_label), np.unique(s_st.results.peak_label))
        st.write('with ' + str(len(s_st.intersection_compounds)) + ' compounds in the intersection' )
        s_st.historical_results = s_st.historical_results
    except:
        st.write('there is a problem with the compounds intersection between the historical data and the current results')
        
        
    ### assessing sample by sample & compound by compound ###
    
    #### TESTING FOR MZ DRIFT ########
    if (('peak_mass_diff_50pc' in s_st.historical_results.columns) == False):
        st.write('#### mz drift cannot be carried out, the columns for mz drift are missing from historical data')
    elif ((('peak_mass_diff_50pc' in s_st.results.columns) == False) & ('peak_mass_diff_50pc' in s_st.historical_results.columns) ):
        st.write('#### mz drift cannot be carried out, the columns for mz are missing from current results data')
    else:
        st.write('running the mz drift analysis')
    try:
        if (('peak_mass_diff_50pc' in s_st.historical_results.columns) & ('peak_mass_diff_50pc' in s_st.results.columns) ):
            st.write('#### indicate a threshold for mz drift')
            s_st.mz_dt = float(st.text_input("maximum acceptable mz drift", '5'))
            
            for compound in s_st.intersection_compounds:
                for sample in s_st.intersection_samples:
                    
                    n1 = np.mean(s_st.historical_results.peak_mass_diff_50pc[(s_st.historical_results.peak_label == compound) & \
                                                                          (s_st.historical_results.STDType == sample)])
                    n2 = np.mean(s_st.results.peak_mass_diff_50pc[(s_st.results.peak_label == compound) & \
                                                                          (s_st.results.STDType == sample)])
                    if abs(n1 - n2) > s_st.mz_dt:
                        st.write('problematic compound: ' + compound + ' in sample: ' + s_st.std_flag +  str(sample) +  ' with ' + \
                                 str(np.round(abs(n1-n2), 2)) + ' ppm drift' )  
                        
    except:
        st.write('there was a problem while running the mz drift analysis')
        
   #### TESTING FOR RT DRIFT ########
    if (('peak_rt_of_max' in s_st.historical_results.columns) == False):
        st.write('#### rt drift cannot be carried out, the columns for rt drift are missing from historical data')
    elif ((('peak_rt_of_max' in s_st.results.columns) == False) & ('peak_rt_of_max' in s_st.historical_results.columns) ):
        st.write('#### mz drift cannot be carried out, the columns for rt drift are missing from current results data')
    else:
        st.write('#### indicate a threshold for RT drift in percent')
        s_st.rt_dt = float( st.text_input("maximum acceptable retention time drift", '1') )
        st.write('#### indicate a threshold for the number of problematic samples')
        s_st.ps_th1 = float( st.text_input("maximum aceptable number of samples for rt drift", '5') )        
        
        st.write('running the rt drift analysis ...')
        
        try:
            if (('peak_rt_of_max' in s_st.historical_results.columns) & ('peak_rt_of_max' in s_st.results.columns) ):
                for compound in s_st.intersection_compounds:
                    st.write(compound)
                    k = 0
                    for sample in s_st.intersection_samples:
                        st.write(sample)
                    
                        n1 = np.mean(s_st.historical_results.peak_rt_of_max[(s_st.historical_results.peak_label == compound) & \
                                                                             (s_st.historical_results.STDType == sample)])
                        n2 = np.mean(s_st.results.peak_rt_of_max[(s_st.results.peak_label == compound) & \
                                                                              (s_st.results.STDType == sample)])

                        st.write(n1)
                        st.write(n2)
                        if abs(n1 - n2) > s_st.rt_dt:
                            k += 1
    #                         st.write('problematic compound: ' + compound + ' in sample: ' +  s_st.std_flag + str(sample) +  ' with ' + \
    #                                  str( np.round(100*abs(n1 - n2)/max(n1,n2), 2)) + ' percent in rt drift' )
                    
                    if k > s_st.ps_th1:
                        st.write(compound + ' showed problems in ' + str(k) + ' samples')
                    else:
                        st.write(compound + ' OK')
                    
        except:
            st.write('there was a problem while running the rt drift analysis')

   #### TESTING FOR PEAK_MAX DRIFT ########
    if (('peak_max' in s_st.historical_results.columns) == False):
        st.write('#### peak intensity drift cannot be carried out, the columns for peak intensity drift are missing from historical data')
    elif ((('peak_max' in s_st.results.columns) == False) & ('peak_max' in s_st.historical_results.columns) ):
        st.write('#### peak intensity drift cannot be carried out, the columns for peak intensity drift are missing from current results data')
    else:
        st.write('#### indicate a threshold for peak height drift in percent')
        s_st.ph_dt = float( st.text_input("peak_max difference in percent", '5') )

        st.write('#### indicate a threshold for the number of problematic samples')
        s_st.ps_th2 = float( st.text_input("maximum aceptable number of samples for peak height drift", '5') )
        
        st.write('running the peak height drift analysis')
    try:
        if (('peak_max' in s_st.historical_results.columns) & ('peak_max' in s_st.results.columns) ):
            
            for compound in s_st.intersection_compounds:
                k = 0
                for sample in s_st.intersection_samples:
                
                    n1 = np.mean(s_st.historical_results.peak_max[(s_st.historical_results.peak_label == compound) & (s_st.historical_results.STDType == sample)])
                    n2 = np.mean(s_st.results.peak_max[(s_st.results.peak_label == compound) & \
                                                                          (s_st.results.STDType == sample)])
                    
                    if abs(n2 - n1)/max(n1,n2) > s_st.ph_dt/100:
                        k += 1
#                         st.write('problematic compound: ' + compound + 'in sample ' + s_st.std_flag + str(sample) + ' with ' +\
#                                  str(100*abs(n2 - n1)/max(n1,n2)) + ' percent in peak height drift' )
                if k > s_st.ps_th2:
                    st.write(compound + ' showed problems in ' + str(k) + ' samples')
                else:
                    st.write(compound + ' OK')
    except:
        st.write('there was a problem while running the peak height drift analysis')
        
st.write(':heavy_minus_sign:' * 35)

st.write(':heavy_minus_sign:' * 35)
st.markdown('''# Growth control section''')



try:

    
    
    st.write('#### indicate the intensity measurement')
    s_st.value_column = st.selectbox('select the intensity measurement \n', ['peak_max', 'peak_area'])
    st.write(s_st.value_column)
    
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
    
    div = np.mean(s_st.results0[s_st.value_column][(s_st.results0.peak_label == s_st.cp) & (s_st.results0.ms_file.str.contains(s_st.gr_flag))]) / \
       np.mean(s_st.results0[s_st.value_column][(s_st.results0.peak_label == s_st.cp) & (s_st.results0.ms_file.str.contains(s_st.me_flag))])
    
    
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
