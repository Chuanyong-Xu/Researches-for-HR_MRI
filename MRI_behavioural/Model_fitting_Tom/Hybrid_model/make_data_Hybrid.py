#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:56:14 2019

@author: pieter
"""

import pandas as pd
import os
import numpy as np

pplist=[1101, 1104, 1105, 1108, 1110, 1114, 1116, 1117, 1120, 1121, 1123, 1124, 1125, 1126, 1132,
        1141, 1142, 1144, 1146, 1147, 1148, 1150, 1155, 1157, 1160, 1161, 1162, 1164, 1165]

column_list = ["Rule", "Stimulus", "Response", "Reward", "Expected value", "PE_estimate"]

for p in pplist:
    os.chdir('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/MRI_behavioural-29subs/Beha_ModCom/data')
    
    filename= "data_29sub_sub{0}.csv".format(int(p))

    data = pd.read_csv(filename,encoding='utf-8')
    
    Rule=data.values[:,1]
    Stim=data.values[:,2]
    #Resp=(data.values[:,5]=='f').astype(int)
    Resp=data.values[:,4]
    Rew=data.values[:,5]
    Rew[Rew==2]=0

    new_filename='Behavioral_data_subject_{0}_Hybrid.csv'.format(int(p))

    new_data=pd.DataFrame({ 'Rule':Rule, 'Stimulus':Stim, 'Response':Resp, 'Reward':Rew, 'Expected value':np.zeros(len(Rule)), 'PE estimate':np.zeros(len(Rule))}, columns = column_list)

    os.chdir('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/MRI_behavioural-29subs/Beha_ModCom/Hybrid_data')
    new_data.to_csv(new_filename, columns = column_list, float_format ='%.3f')
