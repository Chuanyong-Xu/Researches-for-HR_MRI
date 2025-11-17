#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 10:01:43 2019

@author: pieter
"""

import estimation_RW as estimation
import likelihood_RW as likelihood
import sim_data_RW as sim_data

import numpy as np
import pandas as pd
import os

pplist=[1101, 1104, 1105, 1108, 1110, 1114, 1116, 1117, 1120, 1121, 1123, 1124, 1125, 1126, 1132,
        1141, 1142, 1144, 1146, 1147, 1148, 1150, 1155, 1157, 1160, 1161, 1162, 1164, 1165]

method="g"

pars=[]
Lik=[]
loop=0

os.chdir('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/MRI_behavioural-29subs/Beha_ModCom/RW_data')

for p in pplist:
    print("Estimating parameters of subject {}".format(p))
    
    file = 'Behavioral_data_subject_{0}_RW.csv'.format(int(p))
    
    if method=="gradient":
        est=estimation.estimate(file_name=file)
        pars.append(est)
        Lik.append(-est.fun)
    else:
        est = estimation.evol_estimate(file_name=file)
        pars.append(est.x)
        Lik.append(-est.fun)
        print(est.message)
        
    print("Estimated parameters are {} and {}".format(pars[loop][0], pars[loop][1]))
    
    sim_data.simulate_data(alpha = pars[loop][0], beta = pars[loop][1], file_name = file)

    
    print("Estimated the prediction errors" )
    print("LogLikelihood was {}".format(Lik[loop]))
    loop+=1

print(pars)
Learning_rates= [item[0] for item in pars]
Temperatures= [item[1] for item in pars]

print(Learning_rates)
print(Temperatures)
print("Mean learning rate")
print(np.mean(Learning_rates))
print("std learning rate")
print(np.std(Learning_rates))
print("Mean temperature")
print(np.mean(Temperatures))
print("std temperature")
print(np.std(Temperatures))
print("Mean LogLikelihood")
print(np.mean(Lik))
print("std Loglikelihood")
print(np.std(Lik))
AIC=-2*np.asarray(Lik)+ 2*(len(pars[0]))
print("Mean AIC")
print(np.mean(AIC))
print("std AIC")
print(np.std(AIC))

for p in pplist:
    file = 'Behavioral_data_subject_{0}_RW.csv'.format(int(p))
    data = pd.read_csv(file,encoding='utf-8')
    BIC=np.log(len(data.values[:,1]))*len(pars[0])-2*np.asarray(Lik)

new_filename='RW_output.csv'.format(int(p))
column_list = ["Subject", "Learning rate", "Temperature", "LogLik", "AIC", "BIC"]
    
output=pd.DataFrame({ 'Subject':pplist, 'Learning rate':Learning_rates, 'Temperature':Temperatures, 'LogLik':Lik, 'AIC':AIC, 'BIC':BIC}, columns = column_list)
    
output.to_csv(new_filename, columns = column_list, float_format ='%.5f')




