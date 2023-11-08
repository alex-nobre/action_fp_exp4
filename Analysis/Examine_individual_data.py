# -*- coding: utf-8 -*-


import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

# Set seaborn theme to matplotlib
sns.set_theme()

subFileName = input('Enter File Name:\n')
subFile = './Data/Data_files/' + subFileName

# Read data
subData = pd.read_csv(subFile)

# Remove unnecessary columns
subData = subData[['participant', 'date', 'Counterbalance group', 'Handedness', 'block', 'condition', 'orientation', 'foreperiod', 'action_trigger.rt', 'extFixationDuration', 'corrAns', 'Response.keys', 'Response.corr', 'Response.rt']]

# Rename condition columns for clarity
subData=subData.rename(columns={'Response.corr':'Acc'})
subData=subData.rename(columns={'Response.rt':'RT'})
subData=subData.rename(columns={'Counterbalance group':'Counterbalance'})

# Remove empty rows and practice trials
subData = subData[(subData['condition'] != 'practice') & (subData['orientation'].notnull())]

# Create columns for n-1 and n-2 foreperiods by block
subData['oneBackFP'] = subData.groupby(['block'])['foreperiod'].shift(1)

# Check n and % of errors
print('Total errors:')
print(len(subData[subData['Acc'] == 0]))
print((len(subData[subData['Acc'] == 0])/len(subData)) * 100)

# Errors by orientation
print('Errors by orientation:')
print(len(subData[(subData['orientation']=='left') & (subData['Acc']==0)])/
      len(subData[subData['orientation']=='left'])*100)

print(len(subData[(subData['orientation']=='right') & (subData['Acc']==0)])/
      len(subData[subData['orientation']=='right'])*100)
      
# Errors by condition
print('Errors by condition:')
print(len(subData[(subData['condition']=='external') & (subData['Acc']==0)])/
      len(subData[subData['condition']=='external'])*100)

print(len(subData[(subData['condition']=='action') & (subData['Acc']==0)])/
      len(subData[subData['condition']=='action'])*100)
      

# Plot RT by external fixation duration and action trigger latency
plt.figure()
plt.scatter(subData['extFixationDuration'].values, subData['RT'].values, s = 10)
plt.xlabel("External fixation duration (s)")
plt.ylabel("RT (s)")
plt.show()

plt.figure()
plt.scatter(subData['action_trigger.rt'].values, subData['RT'].values, s = 10)
plt.xlabel("Action Trigger Latency (s)")
plt.ylabel("RT (s)")
plt.show()

# Keep only trials with correct responses to analyze RT
subData = subData[(subData['RT'].notnull()) & (subData['Acc'] == 1)]

# Remove outliers
print('notclean: ' + str(len(subData)))
subData = subData[(subData['RT'] < 1) & (subData['RT'] > 0.1)]
print('clean: ' + str(len(subData)))


summaryData=subData.groupby(['foreperiod','condition'],
                           as_index=False)[['RT','Acc']].mean()

summaryPlot=sns.catplot(x="foreperiod", y="RT", kind = 'point',
                          data=summaryData)

plt.show()

plt.figure()
summaryPlot=sns.catplot(x="foreperiod", y="RT", hue = "condition", kind = 'point',
                          data=summaryData)
plt.show()                        
                          
plt.figure()                          
histPlot = sns.histplot(x = "RT", data = subData, bins = 30)                          
                          
plt.show()

