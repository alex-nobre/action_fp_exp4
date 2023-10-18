# -*- coding: utf-8 -*-

# Import libraries
import pandas as pd
import glob
import os

# File paths
filesPath = '.\Data'
gDrivePath = 'G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experiment_4/Analysis/'

# Find files
FileList=glob.glob(filesPath + '/*.txt')
FileList.sort()

nFiles=int(len(FileList))

dataActionFPAll=pd.DataFrame()

for iFile,FileName in enumerate(FileList):
    
    #Get info for this file
    fileInfo = pd.read_table(FileName, nrows = 4, skiprows = [0])
    countBal = fileInfo.iloc[1][1]
    hand = fileInfo.iloc[2][1]
    
    # Read data
    dataActionFP = pd.read_table(FileName, skiprows = [0,1,2,3,4,5,6])
    
    # Remove practice trials
    dataActionFP = dataActionFP[(dataActionFP['FPType'] != 'practice')]
    
    # Create columns for counterbalancing order and handedness
    dataActionFP['countBalance'] = countBal
    dataActionFP['handedness'] = hand
    
    # Create column for FP n-1 by block
    dataActionFP['oneBackFP'] = dataActionFP.groupby(['block'])['FP'].shift(1)
       
    dataActionFPAll=pd.concat([dataActionFPAll, dataActionFP], axis=0)    

# Save to analysis directory
dataActionFPAll.to_csv('./Analysis/'+'dataActionFPAll.csv', index = False)

# Save to google drive for colab use
dataActionFPAll.to_csv(gDrivePath +'dataActionFPAll.csv', index = False)


