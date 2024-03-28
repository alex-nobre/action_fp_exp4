# -*- coding: utf-8 -*-

# Import libraries
import pandas as pd
import glob
import os

# File paths
filesPath = '.\Data\Data_files'
gDrivePath = 'G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experiment_4/Analysis/'

# Find files
FileList=glob.glob(filesPath + '/*.csv')
FileList.sort()

nFiles=int(len(FileList))

dataActionFPAll=pd.DataFrame()

for iFile,FileName in enumerate(FileList):
  
    ID = FileName[18:21]
    
    # Read data
    subData = pd.read_csv(FileName)
    
    # Remove unnecessary columns
    subData = subData[['participant', 'date', 'Counterbalance group', 'Handedness', 'block', 'condition', 'orientation', 'foreperiod', 'action_trigger.rt', 'extFixationDuration', 'ITIDuration', 'corrAns', 'Response.keys', 'Response.corr', 'Response.rt']]
    
    # Rename condition columns for clarity
    subData=subData.rename(columns={'Response.corr':'Acc'})
    subData=subData.rename(columns={'Response.rt':'RT'})
    subData=subData.rename(columns={'Counterbalance group':'Counterbalance'})
    subData=subData.rename(columns={'ITIDuration':'ITI'})
    
    # Remove empty rows and practice trials
    subData = subData[(subData['condition'] != 'practice') & (subData['orientation'].notnull())]

    # Create columns for n-1 and n-2 foreperiods by block
    subData['oneBackFP'] = subData.groupby(['block'])['foreperiod'].shift(1)
    
    # Replace participant's ID by three digits ID from file name
    subData['participant']=ID
    #cols = subData.columns.tolist()
    #cols = cols[-1:] + cols[1:-1]
    #subData = subData[cols]
       
    dataActionFPAll=pd.concat([dataActionFPAll, subData], axis=0)    

# Save to analysis directory
dataActionFPAll.to_csv('./Analysis/'+'dataActionFPAll.csv', index = False)

# Save to google drive for colab use
dataActionFPAll.to_csv(gDrivePath +'dataActionFPAll.csv', index = False)


