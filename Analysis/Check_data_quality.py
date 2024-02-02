# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:33:24 2022

@author: alpno
"""

import pandas as pd
import glob
import os

filesPath = '.\Data\Data_files'

fileList=glob.glob(filesPath + './*.csv')
fileList.sort()

nFiles=int(len(fileList))

qualityTable = pd.DataFrame(columns=['participant','errorRateAll','errorRateLeft','errorRateRight', 
'errorRateAction', 'errorRateExternal', 'trimmedRTRate'])

for file in fileList:
  subData = pd.read_csv(file)
     
  # Rename condition columns for clarity
  subData=subData.rename(columns={'Response.corr':'Acc'})
  subData=subData.rename(columns={'Response.rt':'RT'})
  subData=subData.rename(columns={'Counterbalance group':'Counterbalance'})
     
  # Remove practice trials
  subData = subData[(subData['condition'] != 'practice')]
    
  # ID
  participantID=set(subData['participant'].tolist()).pop()
     
  # Error rate
  errorRateAll=(len(subData[subData['Acc']==0])/len(subData)) * 100 # overall
  errorRateLeft=(len(subData[(subData['Acc']==0)&(subData['orientation']=='left')])/len(subData[subData['orientation']=='left'])) * 100 # left gabors
  errorRateRight=(len(subData[(subData['Acc']==0)&(subData['orientation']=='right')])/len(subData[subData['orientation']=='right'])) * 100 # right gabors
  errorRateAction=(len(subData[(subData['Acc']==0)&(subData['condition']=='action')])/len(subData[subData['condition']=='action'])) * 100 # action condition
  errorRateExternal=(len(subData[(subData['Acc']==0)&(subData['condition']=='external')])/len(subData[subData['condition']=='external'])) * 100 # external condition
     
     
  # Percentage of trials excluded due to RT trimming
  subData=subData[(subData['RT'].notnull())&(subData['Acc']==1)]
  trimmedRTs=subData[(subData['RT']>1)&(subData['RT']<0.15)]     
  trimmedRTRate=(len(trimmedRTs)/len(subData)) * 100
     
  subQualityData=[participantID, errorRateAll, errorRateLeft,
  errorRateRight, errorRateAction, errorRateExternal, 
  trimmedRTRate]
     
  qualityTable.loc[len(qualityTable.index)]=subQualityData
 
qualityTable.to_csv('./Analysis/data_quality.csv')
