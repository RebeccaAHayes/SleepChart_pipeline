# -*- coding: utf-8 -*-
import mne
import yasa
import numpy as np
import pandas as pd
import re
import os
import glob

#set working directory
os.chdir("/lab-share/Psych-Jalbrzikowski-PNR-e2/Public/muse/derivatives/formatted_data")
#os.mkdir("annot_predicted")

#%%
#Define staging function that will return an accuracy in decimal format
def getAccuracy(edf,ch1,fname):
    y=yasa.SleepStaging(edf,eeg_name=ch1)
    y_pred=y.predict()
    y_conf=y.predict_proba()
    y_pred.tofile("annot_predicted/"+fname+"_"+ch1+"_yasa.eannot",sep="\n",format="%s")
    y_conf.to_csv("annot_confidence/"+fname+"_"+ch1+"_yasa.txt",sep="\t")

#%%
#get a list of hypnograms and set up an empty data frame for results
edfs=glob.glob('final/*.edf')

#loop through histograms 
for i in edfs:
    subid=re.search("sub-.*(?=[.])", i).group()
    raw=mne.io.read_raw_edf(i,preload=True)
    getAccuracy(raw,"fp2d_fp1d",subid)
    getAccuracy(raw,"fp1d_fp2d",subid)
    getAccuracy(raw,"fp1d_tp7d",subid)
    getAccuracy(raw,"fp1d_tp8d",subid)
    getAccuracy(raw,"fp2d_tp7d",subid)
    getAccuracy(raw,"fp2d_tp8d",subid)

