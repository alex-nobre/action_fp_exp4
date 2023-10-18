# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:33:24 2022

@author: alpno
"""

import pandas as pd
import glob
import os

filesPath = '.\Data'

fileList=glob.glob(filesPath + './*.txt')
fileList.sort()

nFiles=int(len(fileList))

qualityTable = pd.DataFrame(columns=['participant','errorRateAll','errorRateLeft','errorRateRight', 
'errorRateAction', 'errorRateExternal', 'errorRateConstant', 'errorRateVariable', 'trimmedRTRate'])

for file in fileList:
     subData=pd.read_table(file, skiprows = [0,1,2,3,4,5,6])
     
     # Remove practice trials
     subData = subData[(subData['FPType'] != 'practice')]
    
    # ID
     participantID=set(subData['participant'].tolist()).pop()
     
     # Error rate
     errorRateAll=(len(subData[subData['Acc']==0])/len(subData)) * 100 # overall
     errorRateLeft=(len(subData[(subData['Acc']==0)&(subData['orientation']=='left')])/len(subData[subData['orientation']=='left'])) * 100 # left gabors
     errorRateRight=(len(subData[(subData['Acc']==0)&(subData['orientation']=='right')])/len(subData[subData['orientation']=='right'])) * 100 # right gabors
     errorRateAction=(len(subData[(subData['Acc']==0)&(subData['condition']=='action')])/len(subData[subData['condition']=='action'])) * 100 # action condition
     errorRateExternal=(len(subData[(subData['Acc']==0)&(subData['condition']=='external')])/len(subData[subData['condition']=='external'])) * 100 # external condition
     errorRateConstant=(len(subData[(subData['Acc']==0)&(subData['FPType']=='constant')])/len(subData[subData['FPType']=='constant'])) * 100 # constant FPs
     errorRateVariable=(len(subData[(subData['Acc']==0)&(subData['FPType']=='variable')])/len(subData[subData['FPType']=='variable'])) * 100 # variable FPs
     
     
     # Percentage of trials excluded due to RT trimming
     subData=subData[(subData['RT'].notnull())&(subData['Acc']==1)]
     trimmedRTs=subData[(subData['RT']>1)&(subData['RT']<0.15)]     
     trimmedRTRate=(len(trimmedRTs)/len(subData)) * 100
     
     subQualityData=[participantID, errorRateAll, errorRateLeft,
                     errorRateRight, errorRateAction, errorRateExternal, 
                     errorRateConstant, errorRateVariable, trimmedRTRate]
     
     qualityTable.loc[len(qualityTable.index)]=subQualityData
 
qualityTable.to_csv('./Analysis/data_quality.csv')
