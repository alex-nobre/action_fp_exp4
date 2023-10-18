# -*- coding: utf-8 -*-


import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

# Set seaborn theme to matplotlib
sns.set_theme()

subFileName = input('Enter File Name:\n')
subFile = './Data/' + subFileName

# Read data
subData = pd.read_table(subFile, skiprows = [0,1,2,3,4,5,6])

# Remove practice trials
subData = subData[(subData['FPType'] != 'practice')]

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
      
# Errors by FP Type
print('Errors by FP type:')
print(len(subData[(subData['FPType']=='constant') & (subData['conFPDur']=='short') & (subData['Acc']==0)])/
      len(subData[(subData['FPType']=='constant') & (subData['conFPDur']=='short')])*100)
      
print(len(subData[(subData['FPType']=='constant') & (subData['conFPDur']=='long') & (subData['Acc']==0)])/
      len(subData[(subData['FPType']=='constant') & (subData['conFPDur']=='long')])*100)
      
print(len(subData[(subData['FPType']=='variable') & (subData['Acc']==0)])/
      len(subData[(subData['FPType']=='variable')])*100)

# Inspect n of trials with premature responses
print('N of premature responses:')
len(subData.query('preRespGiven == 1'))

# Plot RT by external fixation duration and action trigger latency
plt.figure()
plt.scatter(subData['extFixDur'].values, subData['RT'].values, s = 10)
plt.show()

plt.figure()
plt.scatter(subData['actionTrigLatency'].values, subData['RT'].values, s = 10)
plt.show()

# Keep only trials with correct responses to analyze RT
subData = subData[(subData['RT'].notnull()) & (subData['Acc'] == 1) & (subData['preRespGiven'] == 0)]

# Remove outliers
print('notclean: ' + str(len(subData)))
subData = subData[(subData['RT'] < 1) & (subData['RT'] > 0.1)]
print('clean: ' + str(len(subData)))


summaryData=subData.groupby(['FP','condition','FPType','block'],
                           as_index=False)[['RT','Acc']].mean()

summaryPlot=sns.catplot(x="FP", y="RT",col = 'FPType', kind = 'point',
                          data=summaryData)
                          
plt.show()
